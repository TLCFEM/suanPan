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

#include "CSMQ.h"

#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>
#include <Toolbox/utility.h>

CSMQ::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , weight(W)
    , m_material(std::move(M)) {}

CSMQ::CSMQ(const unsigned T, uvec&& N, const unsigned M, const unsigned NN, const double TH, const double L)
    : MaterialElement2D(T, NN, m_dof, std::move(N), uvec{M}, false, {DOF::U1, DOF::U2, DOF::UR3})
    , m_node(NN)
    , thickness(TH) { access::rw(characteristic_length) = L; }

int CSMQ::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(!material_proto->is_support_couple()) {
        suanpan_warning("Element {} is assigned with a material that does not support couple stress.\n", get_tag());
        return SUANPAN_FAIL;
    }

    if(PlaneType::E == material_proto->get_plane_type()) suanpan::hacker(thickness) = 1.;

    const auto ele_coor = get_coordinate(2);

    if(characteristic_length < 0.) access::rw(characteristic_length) = sqrt(area::shoelace(ele_coor));

    const IntegrationPlan plan(2, 3, IntegrationType::GAUSS);

    const auto& t_dof = get_translation_dof();
    const auto& r_dof = get_rotation_dof();

    const auto& t_size = t_dof.n_elem;
    const auto& r_size = r_dof.n_elem;
    constexpr auto s_size = 11u;

    mat E1(t_size, t_size, fill::zeros), E2(s_size, s_size, fill::zeros), H1(t_size, t_size, fill::zeros), H2(t_size, s_size, fill::zeros), H3(r_size, t_size, fill::zeros), H4(t_size, t_size, fill::zeros), H5(s_size, s_size, fill::zeros);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1)};
        const auto n = compute_shape_function(t_vec, 0);
        const auto pn = compute_shape_function(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob) * thickness, material_proto->get_copy());

        auto& c_pt = int_pt.back();

        c_pt.m_material->set_characteristic_length(characteristic_length);
        c_pt.m_material->initialize_couple(D);

        mat phi_s(2, t_size, fill::zeros), l_p(3, t_size, fill::zeros), j_p(1, t_size, fill::zeros), j_q(2, r_size, fill::zeros);

        const auto& phi_q = n;
        const auto& phi_r = phi_s;
        const auto& j_s = j_p;

        const mat pnpxy = solve(jacob, pn);
        for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += 2, L += 2) {
            j_q(1, J) = -(j_p(0, L) = l_p(0, K) = l_p(2, L) = pnpxy(0, J));
            j_p(0, K) = -(j_q(0, J) = l_p(2, K) = l_p(1, L) = pnpxy(1, J));
            phi_s(0, K) = phi_s(1, L) = n(J);
        }

        j_p *= .5;
        j_q *= .5;

        const auto location = (n * ele_coor).eval();

        const mat phi_a = shape::stress11(location(0), location(1));
        const mat phi_b = solve(c_pt.m_material->get_initial_stiffness(), phi_a);

        E1 += c_pt.weight * phi_r.t() * c_pt.m_material->get_initial_couple_stiffness() * phi_r;
        E2 += c_pt.weight * phi_b.t() * c_pt.m_material->get_initial_stiffness() * phi_b;
        H1 -= 2. * c_pt.weight * j_p.t() * j_s;
        H2 += c_pt.weight * l_p.t() * phi_a;
        H3 += 2. * c_pt.weight * (phi_q.t() * j_s - j_q.t() * phi_s);
        H4 += c_pt.weight * phi_r.t() * phi_s;
        H5 += c_pt.weight * phi_b.t() * phi_a;

        c_pt.b1 = phi_b;
        c_pt.b2 = c_pt.b3 = phi_r;
    }

    const mat T1 = solve(H5.t(), H2.t());
    const mat T2 = solve(H4.t(), H1.t());
    const mat T3 = solve(H4.t(), H3.t());

    initial_stiffness.set_size(m_size, m_size);
    initial_stiffness(t_dof, t_dof) = T1.t() * E2 * T1 + T2.t() * E1 * T2;
    initial_stiffness(t_dof, r_dof) = T2.t() * E1 * T3;
    initial_stiffness(r_dof, t_dof) = T3.t() * E1 * T2;
    initial_stiffness(r_dof, r_dof) = T3.t() * E1 * T3;

    trial_stiffness = current_stiffness = initial_stiffness;

    for(auto& I : int_pt) {
        I.b1 *= T1;
        I.b2 *= T2;
        I.b3 *= T3;
    }

    if(const auto t_density = material_proto->get_density(); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(const auto& I : int_pt) {
            const auto n_int = compute_shape_function(I.coor, 0);
            const auto t_factor = t_density * I.weight;
            for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof)
                for(auto K = J, M = L; K < m_node; ++K, M += m_dof) initial_mass(L, M) += t_factor * n_int(J) * n_int(K);
        }
        for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    body_force.zeros(m_size, m_dof);
    for(const auto& I : int_pt) {
        const mat n_int = I.weight * compute_shape_function(I.coor, 0);
        for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof)
            for(auto K = 0llu; K < m_dof; ++K) body_force(L + K, K) += n_int(J);
    }

    return SUANPAN_SUCCESS;
}

int CSMQ::update_status() {
    const auto t_disp = get_trial_displacement();

    const auto& t_dof = get_translation_dof();
    const auto& r_dof = get_rotation_dof();

    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    for(const auto& I : int_pt) {
        if(SUANPAN_SUCCESS != I.m_material->update_trial_status(I.b1 * t_disp(t_dof))) return SUANPAN_FAIL;
        if(SUANPAN_SUCCESS != I.m_material->update_couple_trial_status(I.b2 * t_disp(t_dof) + I.b3 * t_disp(r_dof))) return SUANPAN_FAIL;

        trial_stiffness(t_dof, t_dof) += I.weight * I.b1.t() * I.m_material->get_trial_stiffness() * I.b1 + I.weight * I.b2.t() * I.m_material->get_trial_couple_stiffness() * I.b2;
        trial_stiffness(t_dof, r_dof) += I.weight * I.b2.t() * I.m_material->get_trial_couple_stiffness() * I.b3;
        trial_stiffness(r_dof, t_dof) += I.weight * I.b3.t() * I.m_material->get_trial_couple_stiffness() * I.b2;
        trial_stiffness(r_dof, r_dof) += I.weight * I.b3.t() * I.m_material->get_trial_couple_stiffness() * I.b3;

        trial_resistance(t_dof) += I.weight * I.b1.t() * I.m_material->get_trial_stress() + I.weight * I.b2.t() * I.m_material->get_trial_couple_stress();
        trial_resistance(r_dof) += I.weight * I.b3.t() * I.m_material->get_trial_couple_stress();
    }

    return SUANPAN_SUCCESS;
}

int CSMQ::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status() + I.m_material->commit_couple_status();
    return code;
}

int CSMQ::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status() + I.m_material->clear_couple_status();
    return code;
}

int CSMQ::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status() + I.m_material->reset_couple_status();
    return code;
}

mat CSMQ::compute_shape_function(const mat& coordinate, const unsigned order) const { return shape::quad(coordinate, order, m_node); }

std::vector<vec> CSMQ::record(const OutputType P) {
    std::vector<vec> data;
    for(const auto& I : int_pt) append_to(data, I.m_material->record(P));
    return data;
}

void CSMQ::print() {
    suanpan_info("A membrane element (CSMQ) with {} nodes.\n", m_node);
    suanpan_info("The nodes connected are:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("IP {}:\t", I + 1);
        suanpan_info(int_pt[I].coor);
        int_pt[I].m_material->print();
    }
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void CSMQ::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < 4; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void CSMQ::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.head_rows(2) = reshape(get_current_acceleration(), m_dof, m_node).eval().head_rows(2);
    else if(OutputType::V == type) t_disp.head_rows(2) = reshape(get_current_velocity(), m_dof, m_node).eval().head_rows(2);
    else if(OutputType::U == type) t_disp.head_rows(2) = reshape(get_current_displacement(), m_dof, m_node).eval().head_rows(2);

    for(unsigned I = 0; I < 4; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

mat CSMQ::GetData(const OutputType P) {
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
    data.row(4) = interpolation::quadratic(0., -1.);
    data.row(5) = interpolation::quadratic(1., 0.);
    data.row(6) = interpolation::quadratic(0., 1.);
    data.row(7) = interpolation::quadratic(-1., 0.);

    return (data * solve(A, B.t())).t();
}

void CSMQ::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t().eval().head_cols(2);
    for(unsigned I = 0; I < 4; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif

const uvec CSMQ5::t_dof{0, 1, 3, 4, 6, 7, 9, 10, 12, 13};
const uvec CSMQ5::r_dof{2, 5, 8, 11, 14};

CSMQ5::CSMQ5(const unsigned T, uvec&& N, const unsigned M, const double TH, const double L)
    : CSMQ(T, std::move(N), M, 5, TH, L) {}

const uvec& CSMQ5::get_translation_dof() { return t_dof; }

const uvec& CSMQ5::get_rotation_dof() { return r_dof; }

const uvec CSMQ6::t_dof{0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16};
const uvec CSMQ6::r_dof{2, 5, 8, 11, 14, 17};

CSMQ6::CSMQ6(const unsigned T, uvec&& N, const unsigned M, const double TH, const double L)
    : CSMQ(T, std::move(N), M, 6, TH, L) {}

const uvec& CSMQ6::get_translation_dof() { return t_dof; }

const uvec& CSMQ6::get_rotation_dof() { return r_dof; }

const uvec CSMQ7::t_dof{0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19};
const uvec CSMQ7::r_dof{2, 5, 8, 11, 14, 17, 20};

CSMQ7::CSMQ7(const unsigned T, uvec&& N, const unsigned M, const double TH, const double L)
    : CSMQ(T, std::move(N), M, 7, TH, L) {}

const uvec& CSMQ7::get_translation_dof() { return t_dof; }

const uvec& CSMQ7::get_rotation_dof() { return r_dof; }
