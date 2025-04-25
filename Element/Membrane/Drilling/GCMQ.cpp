/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include "GCMQ.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>
#include <Toolbox/tensor.h>
#include <Toolbox/utility.h>

mat GCMQ::form_transformation(const mat& jacobian) {
    mat trans_mat(3, 3);

    trans_mat(0, 0) = jacobian(0, 0) * jacobian(0, 0);
    trans_mat(1, 0) = jacobian(0, 1) * jacobian(0, 1);
    trans_mat(2, 0) = jacobian(0, 0) * jacobian(0, 1);

    trans_mat(0, 1) = jacobian(1, 0) * jacobian(1, 0);
    trans_mat(1, 1) = jacobian(1, 1) * jacobian(1, 1);
    trans_mat(2, 1) = jacobian(1, 0) * jacobian(1, 1);

    trans_mat(0, 2) = 2. * jacobian(0, 0) * jacobian(1, 0);
    trans_mat(1, 2) = 2. * jacobian(1, 0) * jacobian(1, 1);
    trans_mat(2, 2) = jacobian(0, 0) * jacobian(1, 1) + jacobian(0, 1) * jacobian(1, 0);

    return trans_mat / accu(square(trans_mat));
}

mat GCMQ::form_enhanced_strain(const vec& coor, const int num_enhanced_mode) {
    mat poly(3, num_enhanced_mode, fill::zeros);

    if(auto &X = coor(0), &Y = coor(1); 1 == num_enhanced_mode) {
        poly(0, 0) = 3. * X * X - 1.;
        poly(1, 0) = 3. * Y * Y - 1.;
    }
    else if(2 == num_enhanced_mode) {
        poly(2, 1) = poly(0, 0) = 3. * X * X - 1.;
        poly(2, 0) = poly(1, 1) = 3. * Y * Y - 1.;
    }
    else throw invalid_argument("not supported");

    return poly;
}

int GCMQ::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    access::rw(mat_stiffness) = material_proto->get_initial_stiffness();

    if(PlaneType::E == material_proto->get_plane_type()) suanpan::hacker(thickness) = 1.;

    const auto ele_coor = get_coordinate(2);

    access::rw(characteristic_length) = sqrt(area::shoelace(ele_coor));

    access::rw(iso_mapping) = trans(mapping * ele_coor);

    const IntegrationPlan plan(2, scheme == 'I' ? 2 : 3, scheme == 'I' ? IntegrationType::IRONS : scheme == 'L' ? IntegrationType::LOBATTO : IntegrationType::GAUSS);

    const auto diff_coor = form_diff_coor(ele_coor);

    const auto jacob_trans = form_transformation(shape::quad(vec{0., 0.}, 1) * ele_coor);

    mat H(11, 11, fill::zeros), HTT(11, 11, fill::zeros);

    N.zeros(11, 12);
    M.zeros(11, enhanced_mode);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const auto &X = plan(I, 0), &Y = plan(I, 1);

        vec t_vec{X, Y};
        const auto pn = compute_shape_function(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(std::move(t_vec), det(jacob) * plan(I, 2) * thickness, material_proto->get_copy());

        auto& c_pt = int_pt.back();

        c_pt.poly_stress = shape::stress11(iso_mapping * form_stress_mode(X, Y));
        c_pt.poly_strain = solve(mat_stiffness, c_pt.poly_stress);

        M += c_pt.factor * c_pt.poly_stress.t() * jacob_trans * form_enhanced_strain(c_pt.coor, enhanced_mode);
        N += c_pt.factor * c_pt.poly_stress.t() * form_displacement_dn(solve(jacob, pn), solve(jacob, form_drilling_dn(c_pt.coor, diff_coor)));
        H += c_pt.factor * c_pt.poly_stress.t() * c_pt.poly_strain;
        HTT += c_pt.factor * c_pt.poly_strain.t() * mat_stiffness * c_pt.poly_strain;
    }

    HT = trans(H);

    if(!solve(NT, H, N) || !solve(MT, H, M)) {
        suanpan_error("Element {} fails to initialize and is disabled.\n", get_tag());
        return SUANPAN_FAIL;
    }

    const mat T = HTT * MT, W = NT.t() * T;

    trial_stiffness = current_stiffness = initial_stiffness = NT.t() * HTT * NT - W * (trial_viwt = current_viwt = initial_viwt = solve(MT.t() * T, W.t()));

    pre_disp.zeros(m_size);

    trial_vif = current_vif.zeros(enhanced_mode);
    trial_zeta = current_zeta.zeros(enhanced_mode);
    trial_beta = current_beta.zeros(11);
    trial_alpha = current_alpha.zeros(11);
    trial_q = current_q.zeros(11);

    form_mass(material_proto->get_density(), diff_coor);

    form_body_force(diff_coor);

    return SUANPAN_SUCCESS;
}

int GCMQ::update_status() {
    vec incre_disp = -pre_disp;
    incre_disp += pre_disp = get_incre_displacement();

    const vec incre_zeta = -trial_viwt * incre_disp - trial_vif;

    trial_zeta += incre_zeta;
    trial_beta += NT * incre_disp + MT * incre_zeta;

    vec local_stress(11, fill::zeros);
    mat local_stiffness(11, 11, fill::zeros);
    for(const auto& t_pt : int_pt) {
        if(t_pt.m_material->update_trial_status(t_pt.poly_strain * trial_beta) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        local_stress += t_pt.factor * t_pt.poly_strain.t() * t_pt.m_material->get_trial_stress();
        local_stiffness += t_pt.factor * t_pt.poly_strain.t() * t_pt.m_material->get_trial_stiffness() * t_pt.poly_strain;
    }

    const mat T = NT.t() * local_stiffness, V = MT.t() * local_stiffness * MT, W = T * MT;

    trial_alpha = solve(HT, local_stress);

    if(!solve(trial_viwt, V, W.t())) return SUANPAN_FAIL;
    if(!solve(trial_vif, V, M.t() * trial_alpha)) return SUANPAN_FAIL;

    trial_resistance = N.t() * trial_alpha - W * trial_vif;
    trial_stiffness = T * (NT - MT * trial_viwt);

    return SUANPAN_SUCCESS;
}

int GCMQ::commit_status() {
    current_zeta = trial_zeta;
    current_beta = trial_beta;
    current_alpha = trial_alpha;
    current_q = trial_q;
    current_vif = trial_vif;
    current_viwt = trial_viwt;

    pre_disp.zeros();

    return SGCMQ::commit_status();
}

int GCMQ::clear_status() {
    trial_zeta = current_zeta.zeros();
    trial_beta = current_beta.zeros();
    trial_alpha = current_alpha.zeros();
    trial_q = current_q.zeros();
    trial_vif = current_vif.zeros();

    pre_disp.zeros();

    current_viwt = trial_viwt = initial_viwt;

    return SGCMQ::clear_status();
}

int GCMQ::reset_status() {
    trial_zeta = current_zeta;
    trial_beta = current_beta;
    trial_alpha = current_alpha;
    trial_q = current_q;
    trial_vif = current_vif;
    trial_viwt = current_viwt;

    pre_disp.zeros();

    return SGCMQ::reset_status();
}

mat GCMQ::compute_shape_function(const mat& coordinate, const unsigned order) const { return shape::quad(coordinate, order, m_node); }

std::vector<vec> GCMQ::record(const OutputType P) {
    std::vector<vec> data;

    if(P == OutputType::S) for(const auto& I : int_pt) data.emplace_back(I.poly_stress * current_alpha);
    else if(P == OutputType::S11) for(const auto& I : int_pt) data.emplace_back(I.poly_stress.row(0) * current_alpha);
    else if(P == OutputType::S22) for(const auto& I : int_pt) data.emplace_back(I.poly_stress.row(1) * current_alpha);
    else if(P == OutputType::S12) for(const auto& I : int_pt) data.emplace_back(I.poly_stress.row(2) * current_alpha);
    else if(P == OutputType::SP) for(const auto& I : int_pt) data.emplace_back(transform::stress::principal(I.poly_stress * current_alpha));
    else if(P == OutputType::SP1) for(const auto& I : int_pt) data.emplace_back(vec{transform::stress::principal(I.poly_stress * current_alpha).at(0)});
    else if(P == OutputType::SP2) for(const auto& I : int_pt) data.emplace_back(vec{transform::stress::principal(I.poly_stress * current_alpha).at(1)});
    else if(P == OutputType::E) for(const auto& I : int_pt) data.emplace_back(I.poly_strain * current_beta);
    else if(P == OutputType::E11) for(const auto& I : int_pt) data.emplace_back(I.poly_strain.row(0) * current_beta);
    else if(P == OutputType::E22) for(const auto& I : int_pt) data.emplace_back(I.poly_strain.row(1) * current_beta);
    else if(P == OutputType::E12) for(const auto& I : int_pt) data.emplace_back(I.poly_strain.row(2) * current_beta);
    else if(P == OutputType::EP) for(const auto& I : int_pt) data.emplace_back(transform::strain::principal(I.poly_strain * current_beta));
    else if(P == OutputType::EP1) for(const auto& I : int_pt) data.emplace_back(vec{transform::strain::principal(I.poly_strain * current_beta).at(0)});
    else if(P == OutputType::EP2) for(const auto& I : int_pt) data.emplace_back(vec{transform::strain::principal(I.poly_strain * current_beta).at(1)});
    else if(P == OutputType::PE) for(const auto& I : int_pt) { data.emplace_back(I.poly_strain * current_beta - solve(mat_stiffness, I.poly_stress * current_alpha)); }
    else if(P == OutputType::PEP) for(const auto& I : int_pt) data.emplace_back(transform::strain::principal(I.poly_strain * current_beta - solve(mat_stiffness, I.poly_stress * current_alpha)));
    else for(const auto& I : int_pt) append_to(data, I.m_material->record(P));

    return data;
}

void GCMQ::print() {
    suanpan_info("A GCMQ mixed quadrilateral element connecting nodes:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material Response:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("IP {}:\t", I + 1);
        suanpan_info(int_pt[I].coor);
        int_pt[I].m_material->print();
    }
    suanpan_info("Element Response:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("IP {}:\t", I + 1);
        suanpan_info(int_pt[I].coor);
        suanpan_info("Strain:\t", vec{int_pt[I].poly_strain * current_beta});
        suanpan_info("Stress:\t", vec{int_pt[I].poly_stress * current_alpha});
    }
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

mat GCMQ::GetData(const OutputType P) {
    if(OutputType::S == P) {
        mat t_stress(6, m_node, fill::zeros);
        t_stress(uvec{0, 1, 3}, uvec{0}) = shape::stress11(iso_mapping * form_stress_mode(-1., -1.)) * current_alpha;
        t_stress(uvec{0, 1, 3}, uvec{1}) = shape::stress11(iso_mapping * form_stress_mode(1., -1.)) * current_alpha;
        t_stress(uvec{0, 1, 3}, uvec{2}) = shape::stress11(iso_mapping * form_stress_mode(1., 1.)) * current_alpha;
        t_stress(uvec{0, 1, 3}, uvec{3}) = shape::stress11(iso_mapping * form_stress_mode(-1., 1.)) * current_alpha;
        return t_stress;
    }
    if(OutputType::E == P) {
        mat t_strain(6, m_node, fill::zeros);
        t_strain(uvec{0, 1, 3}, uvec{0}) = solve(mat_stiffness, shape::stress11(iso_mapping * form_stress_mode(-1., -1.)) * current_beta);
        t_strain(uvec{0, 1, 3}, uvec{1}) = solve(mat_stiffness, shape::stress11(iso_mapping * form_stress_mode(1., -1.)) * current_beta);
        t_strain(uvec{0, 1, 3}, uvec{2}) = solve(mat_stiffness, shape::stress11(iso_mapping * form_stress_mode(1., 1.)) * current_beta);
        t_strain(uvec{0, 1, 3}, uvec{3}) = solve(mat_stiffness, shape::stress11(iso_mapping * form_stress_mode(-1., 1.)) * current_beta);
        return t_strain;
    }

    mat A(int_pt.size(), 9);
    mat B(6, int_pt.size(), fill::zeros);

    for(size_t I = 0; I < int_pt.size(); ++I) {
        if(const auto C = int_pt[I].m_material->record(P); !C.empty()) B(0, I, size(C[0])) = C[0];
        A.row(I) = interpolation::quadratic(int_pt[I].coor);
    }

    mat data(m_node, 9);

    data.row(0) = interpolation::quadratic(-1., -1.);
    data.row(1) = interpolation::quadratic(1., -1.);
    data.row(2) = interpolation::quadratic(1., 1.);
    data.row(3) = interpolation::quadratic(-1., 1.);

    return (data * solve(A, B.t())).t();
}

#endif
