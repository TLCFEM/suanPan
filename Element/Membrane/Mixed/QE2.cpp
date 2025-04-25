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

#include "QE2.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>
#include <Toolbox/utility.h>

mat QE2::mapping = [] {
    mat t_mapping(4, 4);
    t_mapping.fill(.25);
    t_mapping(1, 0) = t_mapping(1, 3) = t_mapping(2, 0) = t_mapping(2, 1) = t_mapping(3, 1) = t_mapping(3, 3) = -.25;
    return t_mapping;
}();

QE2::IntegrationPoint::IntegrationPoint(vec&& C, const double F, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , factor(F)
    , m_material(std::move(M))
    , B(3, m_size, fill::zeros)
    , BI(3, 2, fill::zeros) {}

vec QE2::form_stress_mode(const double X, const double Y) { return vec{0., X, Y, X * Y}; }

QE2::QE2(const unsigned T, uvec&& N, const unsigned M, const double TH)
    : MaterialElement2D(T, m_node, m_dof, std::move(N), uvec{M}, false, {DOF::U1, DOF::U2})
    , thickness(TH) {}

int QE2::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get<Material>(material_tag(0));

    if(PlaneType::E == mat_proto->get_plane_type()) suanpan::hacker(thickness) = 1.;

    access::rw(mat_stiffness) = mat_proto->get_initial_stiffness();

    const auto ele_coor = get_coordinate(2);

    access::rw(characteristic_length) = sqrt(area::shoelace(ele_coor));

    access::rw(iso_mapping) = trans(mapping * ele_coor);

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

    mat H(7, 7, fill::zeros);

    L.zeros(7, 8);
    LI.zeros(7, 2);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1)};
        const auto pn = compute_shape_function(t_vec, 1);
        const mat jacob = pn * ele_coor;
        const auto det_jacob = det(jacob);
        int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det_jacob * thickness, mat_proto->get_copy());

        auto& c_pt = int_pt.back();

        const auto &X = c_pt.coor(0), &Y = c_pt.coor(1);

        c_pt.A = solve(mat_stiffness, c_pt.P = shape::stress7(iso_mapping * form_stress_mode(X, Y)));

        const mat pn_pxy = solve(jacob, pn);
        for(auto J = 0u, K = 0u, M = 1u; J < m_node; ++J, K += m_dof, M += m_dof) {
            c_pt.B(0, K) = c_pt.B(2, M) = pn_pxy(0, J);
            c_pt.B(2, K) = c_pt.B(1, M) = pn_pxy(1, J);
        }

        c_pt.BI(2, 1) = c_pt.BI(0, 0) = (-iso_mapping(1, 2) * X - iso_mapping(1, 1) * Y) / det_jacob;
        c_pt.BI(2, 0) = c_pt.BI(1, 1) = (iso_mapping(0, 2) * X + iso_mapping(0, 1) * Y) / det_jacob;

        H += c_pt.P.t() * c_pt.A * c_pt.factor;
        L += c_pt.P.t() * c_pt.B * c_pt.factor;
        LI += c_pt.P.t() * c_pt.BI * c_pt.factor;
    }

    HT = trans(H);

    if(!solve(HIL, H, L) || !solve(HILI, H, LI)) {
        suanpan_error("Element {} fails to initialize and is disabled.\n", get_tag());
        return SUANPAN_FAIL;
    }

    const mat TT = LI.t() * HIL;

    trial_stiffness = current_stiffness = initial_stiffness = HIL.t() * L - TT.t() * (trial_qtitt = current_qtitt = initial_qtitt = solve(mat(LI.t() * HILI), TT));

    pre_disp.zeros(m_size);

    current_qtifi = trial_qtifi.zeros(2);
    current_lambda = trial_lambda.zeros(2);
    current_alpha = trial_alpha.zeros(7);
    current_beta = trial_beta.zeros(7);

    if(const auto t_density = mat_proto->get_density(); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(const auto& I : int_pt) {
            const auto n_int = compute_shape_function(I.coor, 0);
            const auto t_factor = t_density * I.factor;
            for(auto J = 0u, P = 0u; J < m_node; ++J, P += m_dof) for(auto K = J, M = P; K < m_node; ++K, M += m_dof) initial_mass(P, M) += t_factor * n_int(J) * n_int(K);
        }
        for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, P = K + m_dof; J < m_size; J += m_dof, P += m_dof) initial_mass(J, I) = initial_mass(K, P) = initial_mass(P, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    body_force.zeros(m_size, m_dof);
    for(const auto& I : int_pt) {
        const mat n_int = I.factor * compute_shape_function(I.coor, 0);
        for(auto J = 0u, M = 0u; J < m_node; ++J, M += m_dof) for(auto K = 0llu; K < m_dof; ++K) body_force(M + K, K) += n_int(J);
    }

    return SUANPAN_SUCCESS;
}

int QE2::update_status() {
    vec incre_disp = -pre_disp;
    incre_disp += pre_disp = get_incre_displacement();

    const vec incre_lambda = -trial_qtitt * incre_disp - trial_qtifi; // eq. 65

    trial_lambda += incre_lambda;                          // eq. 66
    trial_alpha += HIL * incre_disp + HILI * incre_lambda; // eq. 46

    vec local_stress(7, fill::zeros);
    mat local_stiffness(7, 7, fill::zeros);
    for(const auto& t_pt : int_pt) {
        if(t_pt.m_material->update_trial_status(t_pt.A * trial_alpha) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        local_stiffness += t_pt.factor * t_pt.A.t() * t_pt.m_material->get_trial_stiffness() * t_pt.A; // eq. 56
        local_stress += t_pt.factor * t_pt.A.t() * t_pt.m_material->get_trial_stress();
    }

    const mat HILIHT = HILI.t() * local_stiffness, QT = HILIHT * HILI, TT = HILIHT * HIL; // eq. 60

    if(!solve(trial_beta, HT, local_stress)) return SUANPAN_FAIL;
    if(!solve(trial_qtitt, QT, TT)) return SUANPAN_FAIL;                  // eq. 65
    if(!solve(trial_qtifi, QT, LI.t() * trial_beta)) return SUANPAN_FAIL; // eq. 65

    trial_resistance = L.t() * trial_beta - TT.t() * trial_qtifi;             // eq. 64
    trial_stiffness = HIL.t() * local_stiffness * HIL - TT.t() * trial_qtitt; // eq. 61

    return SUANPAN_SUCCESS;
}

int QE2::clear_status() {
    current_lambda = trial_lambda.zeros();
    current_alpha = trial_alpha.zeros();
    current_beta = trial_beta.zeros();
    current_qtifi = trial_qtifi.zeros();
    pre_disp.zeros();

    current_qtitt = trial_qtitt = initial_qtitt;

    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int QE2::commit_status() {
    current_lambda = trial_lambda;
    current_alpha = trial_alpha;
    current_beta = trial_beta;
    current_qtifi = trial_qtifi;
    current_qtitt = trial_qtitt;
    pre_disp.zeros();

    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int QE2::reset_status() {
    trial_lambda = current_lambda;
    trial_alpha = current_alpha;
    trial_beta = current_beta;
    trial_qtifi = current_qtifi;
    trial_qtitt = current_qtitt;
    pre_disp.zeros();

    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

mat QE2::compute_shape_function(const mat& coordinate, const unsigned order) const { return shape::quad(coordinate, order, m_node); }

std::vector<vec> QE2::record(const OutputType P) {
    std::vector<vec> data;

    if(P == OutputType::E) for(const auto& I : int_pt) data.emplace_back(I.A * current_alpha);
    else if(P == OutputType::S) for(const auto& I : int_pt) data.emplace_back(I.P * current_beta);
    else for(const auto& I : int_pt) append_to(data, I.m_material->record(P));

    return data;
}

void QE2::print() {
    suanpan_info("Piltner's mixed quad element connects nodes:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material Response:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("IP {}:\t", I + 1);
        int_pt[I].m_material->print();
    }
    suanpan_info("Element Response:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("IP {}:\n", I + 1);
        suanpan_info("Strain:\t", vec{int_pt[I].A * current_alpha});
        suanpan_info("Stress:\t", vec{int_pt[I].P * current_beta});
    }
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void QE2::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void QE2::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration(), m_dof, m_node);
    else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity(), m_dof, m_node);
    else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

mat QE2::GetData(const OutputType P) {
    if(OutputType::S == P) {
        mat t_stress(6, m_node, fill::zeros);
        t_stress(uvec{0, 1, 3}, uvec{0}) = shape::stress7(iso_mapping * form_stress_mode(-1., -1.)) * current_alpha;
        t_stress(uvec{0, 1, 3}, uvec{1}) = shape::stress7(iso_mapping * form_stress_mode(1., -1.)) * current_alpha;
        t_stress(uvec{0, 1, 3}, uvec{2}) = shape::stress7(iso_mapping * form_stress_mode(1., 1.)) * current_alpha;
        t_stress(uvec{0, 1, 3}, uvec{3}) = shape::stress7(iso_mapping * form_stress_mode(-1., 1.)) * current_alpha;
        return t_stress;
    }
    if(OutputType::E == P) {
        mat t_stress(6, m_node, fill::zeros);
        t_stress(uvec{0, 1, 3}, uvec{0}) = solve(mat_stiffness, shape::stress7(iso_mapping * form_stress_mode(-1., -1.)) * current_beta);
        t_stress(uvec{0, 1, 3}, uvec{1}) = solve(mat_stiffness, shape::stress7(iso_mapping * form_stress_mode(1., -1.)) * current_beta);
        t_stress(uvec{0, 1, 3}, uvec{2}) = solve(mat_stiffness, shape::stress7(iso_mapping * form_stress_mode(1., 1.)) * current_beta);
        t_stress(uvec{0, 1, 3}, uvec{3}) = solve(mat_stiffness, shape::stress7(iso_mapping * form_stress_mode(-1., 1.)) * current_beta);
        return t_stress;
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

void QE2::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
