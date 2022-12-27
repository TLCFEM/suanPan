/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include "CSMT6.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/utility.h>

CSMT6::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::forward<vec>(C))
    , weight(W)
    , m_material(std::forward<unique_ptr<Material>>(M)) {}

CSMT6::CSMT6(const unsigned T, uvec&& N, const unsigned M, const double TH, const double L)
    : MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, false, {DOF::U1, DOF::U2, DOF::UR3})
    , thickness(TH) { access::rw(characteristic_length) = L; }

int CSMT6::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(!material_proto->is_support_couple()) {
        suanpan_warning("Element %u is assigned with a material that does not support couple stress.\n", get_tag());
        return SUANPAN_FAIL;
    }

    if(PlaneType::E == static_cast<PlaneType>(material_proto->get_parameter(ParameterType::PLANETYPE))) suanpan::hacker(thickness) = 1.;

    mat ele_coor(m_node, m_node, fill::none);
    ele_coor.col(0).fill(1.);
    ele_coor.cols(1, 2) = get_coordinate(2);
    ele_coor.col(3) = ele_coor.col(1) % ele_coor.col(2);
    ele_coor.col(4) = square(ele_coor.col(1));
    ele_coor.col(5) = square(ele_coor.col(2));

    const mat inv_coor = inv(ele_coor);

    access::rw(area) = .5 * det(ele_coor(span(0, 2), span(0, 2)));

    if(characteristic_length < 0.) access::rw(characteristic_length) = sqrt(area);

    const IntegrationPlan plan(2, 5, IntegrationType::TRIANGLE);

    const auto& t_size = t_dof.n_elem;
    const auto& r_size = r_dof.n_elem;

    mat E1(t_size, t_size, fill::zeros), E2(9, 9, fill::zeros), H1(t_size, t_size, fill::zeros), H2(t_size, 9, fill::zeros), H3(r_size, t_size, fill::zeros), H4(t_size, t_size, fill::zeros), H5(9, 9, fill::zeros);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec = ele_coor(span(0, 2), span(1, 2)).t() * vec{plan(I, 0), plan(I, 1), plan(I, 2)};
        const rowvec n = shape::triangle(t_vec, 0) * inv_coor;
        const mat pnpxy = shape::triangle(t_vec, 1) * inv_coor;
        int_pt.emplace_back(std::move(t_vec), area * plan(I, 3) * thickness, material_proto->get_copy());

        auto& c_pt = int_pt.back();

        c_pt.m_material->set_characteristic_length(characteristic_length);
        c_pt.m_material->initialize_couple(D);

        mat phi_s(2, t_size, fill::zeros), l_p(3, t_size, fill::zeros), j_p(1, t_size, fill::zeros), j_q(2, r_size, fill::zeros);

        const auto& phi_q = n;
        const auto& phi_r = phi_s;
        const auto& j_s = j_p;

        for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += 2, L += 2) {
            j_q(1, J) = -(j_p(0, L) = l_p(0, K) = l_p(2, L) = pnpxy(0, J));
            j_p(0, K) = -(j_q(0, J) = l_p(2, K) = l_p(1, L) = pnpxy(1, J));
            phi_s(0, K) = phi_s(1, L) = n(J);
        }

        j_p *= .5;
        j_q *= .5;

        const mat location = (n * ele_coor.cols(1, 2)).eval();

        const mat phi_a = shape::linear_stress(location(0), location(1));
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

    if(const auto t_density = material_proto->get_parameter(ParameterType::DENSITY); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(const auto& I : int_pt) {
            const mat n_int = shape::triangle(I.coor, 0) * inv_coor;
            const auto t_factor = t_density * I.weight;
            for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = J, M = L; K < m_node; ++K, M += m_dof) initial_mass(L, M) += t_factor * n_int(J) * n_int(K);
        }
        for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    body_force.zeros(m_size, m_dof);
    for(const auto& I : int_pt) {
        const mat n_int = I.weight * shape::triangle(I.coor, 0) * inv_coor;
        for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = 0llu; K < m_dof; ++K) body_force(L + K, K) += n_int(J);
    }

    return SUANPAN_SUCCESS;
}

int CSMT6::update_status() {
    const auto t_disp = get_trial_displacement();

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

int CSMT6::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status() + I.m_material->commit_couple_status();
    return code;
}

int CSMT6::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status() + I.m_material->clear_couple_status();
    return code;
}

int CSMT6::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status() + I.m_material->reset_couple_status();
    return code;
}

vector<vec> CSMT6::record(const OutputType P) {
    vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.m_material->record(P)) data.emplace_back(J);
    return data;
}

void CSMT6::print() {
    suanpan_info("Element %u is a six-node triangular membrane element (CSMT6).\n", get_tag());
    node_encoding.t().print("The nodes connected are:");
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("Integration Point %llu:\t", I + 1);
        int_pt[I].coor.t().print();
        int_pt[I].m_material->print();
    }
}

#ifdef SUANPAN_VTK
#include <vtkQuadraticTriangle.h>

void CSMT6::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuadraticTriangle>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void CSMT6::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_data(6, m_node, fill::zeros);

    if(OutputType::A == type) t_data.rows(0, 1) = reshape(get_current_acceleration(), m_dof, m_node).eval().head_rows(2);
    else if(OutputType::V == type) t_data.rows(0, 1) = reshape(get_current_velocity(), m_dof, m_node).eval().head_rows(2);
    else if(OutputType::U == type) t_data.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node).eval().head_rows(2);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_data.colptr(I));
}

void CSMT6::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).eval().head_rows(2).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
