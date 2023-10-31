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

#include "PS.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>
#include <Toolbox/utility.h>

PS::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::forward<vec>(C))
    , weight(W)
    , m_material(std::forward<unique_ptr<Material>>(M)) {}

mat PS::form_transformation(const mat& jacobian) {
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

    return trans_mat;
}

PS::PS(const unsigned T, uvec&& N, const unsigned M, const double TH)
    : MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, false, {DOF::U1, DOF::U2})
    , thickness(TH) {}

int PS::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(PlaneType::E == material_proto->get_plane_type()) suanpan::hacker(thickness) = 1.;

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const auto ele_coor = get_coordinate(2);

    access::rw(characteristic_length) = sqrt(area::shoelace(ele_coor));

    const auto jacob_trans = form_transformation(shape::quad(vec{0., 0.}, 1) * ele_coor);

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

    mat poly_disp(3, 8, fill::zeros), poly_stress(3, 5, fill::zeros);
    for(auto J = 0; J < 3; ++J) poly_stress(J, J) = 1.;

    mat H(5, 5, fill::zeros), N(5, 8, fill::zeros);
    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1)};
        const auto pn = compute_shape_function(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(std::move(t_vec), thickness * plan(I, 2) * det(jacob), material_proto->get_copy());

        auto& c_pt = int_pt.back();

        const mat pn_pxy = solve(jacob, pn);
        for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
            poly_disp(2, L) = poly_disp(0, K) = pn_pxy(0, J);
            poly_disp(2, K) = poly_disp(1, L) = pn_pxy(1, J);
        }

        poly_stress.col(3) = jacob_trans.col(0) * c_pt.coor(1);
        poly_stress.col(4) = jacob_trans.col(1) * c_pt.coor(0);

        c_pt.poly_strain = solve(ini_stiffness, poly_stress);

        N += c_pt.weight * poly_stress.t() * poly_disp;
        H += c_pt.weight * poly_stress.t() * c_pt.poly_strain;
    }

    const mat NT = solve(H, N);

    trial_stiffness = current_stiffness = initial_stiffness = N.t() * NT;

    for(auto& I : int_pt) I.poly_strain *= NT;

    if(const auto t_density = material_proto->get_density(); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(const auto& I : int_pt) {
            const auto n_int = compute_shape_function(I.coor, 0);
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
        const mat n_int = I.weight * compute_shape_function(I.coor, 0);
        for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = 0llu; K < m_dof; ++K) body_force(L + K, K) += n_int(J);
    }

    return SUANPAN_SUCCESS;
}

int PS::update_status() {
    const auto trial_disp = get_trial_displacement();

    trial_resistance.zeros(m_size);
    trial_stiffness.zeros(m_size, m_size);
    for(const auto& I : int_pt) {
        if(I.m_material->update_trial_status(I.poly_strain * trial_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        trial_resistance += I.weight * I.poly_strain.t() * I.m_material->get_trial_stress();
        trial_stiffness += I.weight * I.poly_strain.t() * I.m_material->get_trial_stiffness() * I.poly_strain;
    }

    return SUANPAN_SUCCESS;
}

int PS::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int PS::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int PS::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

mat PS::compute_shape_function(const mat& coordinate, const unsigned order) const { return shape::quad(coordinate, order, m_node); }

vector<vec> PS::record(const OutputType P) {
    vector<vec> data;
    for(const auto& I : int_pt) append_to(data, I.m_material->record(P));
    return data;
}

void PS::print() {
    suanpan_info("A four-node membrane element (Pian-Sumihara) connecting nodes:", node_encoding);
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

void PS::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void PS::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration(), m_dof, m_node);
    else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity(), m_dof, m_node);
    else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void PS::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
