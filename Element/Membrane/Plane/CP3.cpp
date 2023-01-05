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

#include "CP3.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/tensor.h>
#include <Toolbox/utility.h>

void CP3::stack_stiffness(mat& K, const mat& D, const mat& N, const double F) {
    const auto D11 = F * D(0, 0);
    const auto D12 = F * D(0, 1);
    const auto D13 = F * D(0, 2);
    const auto D21 = F * D(1, 0);
    const auto D22 = F * D(1, 1);
    const auto D23 = F * D(1, 2);
    const auto D31 = F * D(2, 0);
    const auto D32 = F * D(2, 1);
    const auto D33 = F * D(2, 2);

    const auto& NX1 = N(0, 0);
    const auto& NY1 = N(1, 1);
    const auto& NX2 = N(0, 2);
    const auto& NY2 = N(1, 3);
    const auto& NX3 = N(0, 4);
    const auto& NY3 = N(1, 5);

    const auto D11NX1 = D11 * NX1;
    const auto D11NX2 = D11 * NX2;
    const auto D11NX3 = D11 * NX3;

    const auto D12NX1 = D12 * NX1;
    const auto D12NX2 = D12 * NX2;
    const auto D12NX3 = D12 * NX3;

    const auto D13NX1 = D13 * NX1;
    const auto D13NX2 = D13 * NX2;
    const auto D13NX3 = D13 * NX3;

    const auto D21NY1 = D21 * NY1;
    const auto D21NY2 = D21 * NY2;
    const auto D21NY3 = D21 * NY3;

    const auto D22NY1 = D22 * NY1;
    const auto D22NY2 = D22 * NY2;
    const auto D22NY3 = D22 * NY3;

    const auto D23NY1 = D23 * NY1;
    const auto D23NY2 = D23 * NY2;
    const auto D23NY3 = D23 * NY3;

    const auto D31NX1 = D31 * NX1;
    const auto D31NX2 = D31 * NX2;
    const auto D31NX3 = D31 * NX3;
    const auto D31NY1 = D31 * NY1;
    const auto D31NY2 = D31 * NY2;
    const auto D31NY3 = D31 * NY3;

    const auto D32NX1 = D32 * NX1;
    const auto D32NX2 = D32 * NX2;
    const auto D32NX3 = D32 * NX3;
    const auto D32NY1 = D32 * NY1;
    const auto D32NY2 = D32 * NY2;
    const auto D32NY3 = D32 * NY3;

    const auto D33NX1 = D33 * NX1;
    const auto D33NX2 = D33 * NX2;
    const auto D33NX3 = D33 * NX3;
    const auto D33NY1 = D33 * NY1;
    const auto D33NY2 = D33 * NY2;
    const auto D33NY3 = D33 * NY3;

    const auto D11NX1D31NY1 = D11NX1 + D31NY1;
    const auto D13NX1D33NY1 = D13NX1 + D33NY1;
    const auto D12NX1D32NY1 = D12NX1 + D32NY1;
    const auto D31NX1D21NY1 = D31NX1 + D21NY1;
    const auto D33NX1D23NY1 = D33NX1 + D23NY1;
    const auto D32NX1D22NY1 = D32NX1 + D22NY1;
    const auto D11NX2D31NY2 = D11NX2 + D31NY2;
    const auto D13NX2D33NY2 = D13NX2 + D33NY2;
    const auto D12NX2D32NY2 = D12NX2 + D32NY2;
    const auto D31NX2D21NY2 = D31NX2 + D21NY2;
    const auto D33NX2D23NY2 = D33NX2 + D23NY2;
    const auto D32NX2D22NY2 = D32NX2 + D22NY2;
    const auto D11NX3D31NY3 = D11NX3 + D31NY3;
    const auto D13NX3D33NY3 = D13NX3 + D33NY3;
    const auto D12NX3D32NY3 = D12NX3 + D32NY3;
    const auto D31NX3D21NY3 = D31NX3 + D21NY3;
    const auto D33NX3D23NY3 = D33NX3 + D23NY3;
    const auto D32NX3D22NY3 = D32NX3 + D22NY3;

    K.set_size(m_size, m_size);
    K(0, 0) = NX1 * D11NX1D31NY1 + NY1 * D13NX1D33NY1;
    K(0, 1) = NX1 * D13NX1D33NY1 + NY1 * D12NX1D32NY1;
    K(0, 2) = NX2 * D11NX1D31NY1 + NY2 * D13NX1D33NY1;
    K(0, 3) = NX2 * D13NX1D33NY1 + NY2 * D12NX1D32NY1;
    K(0, 4) = NX3 * D11NX1D31NY1 + NY3 * D13NX1D33NY1;
    K(0, 5) = NX3 * D13NX1D33NY1 + NY3 * D12NX1D32NY1;
    K(1, 0) = NX1 * D31NX1D21NY1 + NY1 * D33NX1D23NY1;
    K(1, 1) = NX1 * D33NX1D23NY1 + NY1 * D32NX1D22NY1;
    K(1, 2) = NX2 * D31NX1D21NY1 + NY2 * D33NX1D23NY1;
    K(1, 3) = NX2 * D33NX1D23NY1 + NY2 * D32NX1D22NY1;
    K(1, 4) = NX3 * D31NX1D21NY1 + NY3 * D33NX1D23NY1;
    K(1, 5) = NX3 * D33NX1D23NY1 + NY3 * D32NX1D22NY1;
    K(2, 0) = NX1 * D11NX2D31NY2 + NY1 * D13NX2D33NY2;
    K(2, 1) = NX1 * D13NX2D33NY2 + NY1 * D12NX2D32NY2;
    K(2, 2) = NX2 * D11NX2D31NY2 + NY2 * D13NX2D33NY2;
    K(2, 3) = NX2 * D13NX2D33NY2 + NY2 * D12NX2D32NY2;
    K(2, 4) = NX3 * D11NX2D31NY2 + NY3 * D13NX2D33NY2;
    K(2, 5) = NX3 * D13NX2D33NY2 + NY3 * D12NX2D32NY2;
    K(3, 0) = NX1 * D31NX2D21NY2 + NY1 * D33NX2D23NY2;
    K(3, 1) = NX1 * D33NX2D23NY2 + NY1 * D32NX2D22NY2;
    K(3, 2) = NX2 * D31NX2D21NY2 + NY2 * D33NX2D23NY2;
    K(3, 3) = NX2 * D33NX2D23NY2 + NY2 * D32NX2D22NY2;
    K(3, 4) = NX3 * D31NX2D21NY2 + NY3 * D33NX2D23NY2;
    K(3, 5) = NX3 * D33NX2D23NY2 + NY3 * D32NX2D22NY2;
    K(4, 0) = NX1 * D11NX3D31NY3 + NY1 * D13NX3D33NY3;
    K(4, 1) = NX1 * D13NX3D33NY3 + NY1 * D12NX3D32NY3;
    K(4, 2) = NX2 * D11NX3D31NY3 + NY2 * D13NX3D33NY3;
    K(4, 3) = NX2 * D13NX3D33NY3 + NY2 * D12NX3D32NY3;
    K(4, 4) = NX3 * D11NX3D31NY3 + NY3 * D13NX3D33NY3;
    K(4, 5) = NX3 * D13NX3D33NY3 + NY3 * D12NX3D32NY3;
    K(5, 0) = NX1 * D31NX3D21NY3 + NY1 * D33NX3D23NY3;
    K(5, 1) = NX1 * D33NX3D23NY3 + NY1 * D32NX3D22NY3;
    K(5, 2) = NX2 * D31NX3D21NY3 + NY2 * D33NX3D23NY3;
    K(5, 3) = NX2 * D33NX3D23NY3 + NY2 * D32NX3D22NY3;
    K(5, 4) = NX3 * D31NX3D21NY3 + NY3 * D33NX3D23NY3;
    K(5, 5) = NX3 * D33NX3D23NY3 + NY3 * D32NX3D22NY3;
}

CP3::CP3(const unsigned T, uvec&& NT, const unsigned MT, const double TH, const bool R)
    : MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(NT), uvec{MT}, R, {DOF::U1, DOF::U2})
    , thickness(TH) {}

int CP3::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(PlaneType::E == static_cast<PlaneType>(material_proto->get_parameter(ParameterType::PLANETYPE))) suanpan::hacker(thickness) = 1.;

    m_material = material_proto->get_copy();

    mat ele_coor(m_node, m_node);
    ele_coor.col(0).fill(1.);
    ele_coor.cols(1, 2) = get_coordinate(2);

    access::rw(area) = .5 * det(ele_coor);

    access::rw(characteristic_length) = sqrt(area);

    const mat inv_coor = inv(ele_coor);
    pn_pxy = inv_coor.rows(1, 2);

    strain_mat.zeros(3, m_size);
    for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
        strain_mat(0, K) = strain_mat(2, L) = pn_pxy(0, J);
        strain_mat(2, K) = strain_mat(1, L) = pn_pxy(1, J);
    }
    trial_stiffness = current_stiffness = initial_stiffness = area * thickness * strain_mat.t() * m_material->get_initial_stiffness() * strain_mat;

    rowvec n = mean(ele_coor) * inv_coor;

    if(const auto t_density = area * thickness * m_material->get_parameter(ParameterType::DENSITY); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(auto I = 0u, K = 0u; I < m_node; ++I, K += m_dof) for(auto J = I, L = K; J < m_node; ++J, L += m_dof) initial_mass(K, L) += t_density * n(I) * n(J);
        for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    body_force.zeros(m_size, m_dof);
    n *= area * thickness;
    for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = 0llu; K < m_dof; ++K) body_force(L + K, K) = n(J);

    return SUANPAN_SUCCESS;
}

int CP3::update_status() {
    if(nlgeom) {
        trial_geometry.zeros(m_size, m_size);

        const mat gradient = reshape(get_trial_displacement(), m_dof, m_node) * pn_pxy.t() + eye(m_dof, m_dof);

        mat BN(3, m_size);
        for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
            BN(0, K) = pn_pxy(0, J) * gradient(0, 0);
            BN(1, K) = pn_pxy(1, J) * gradient(0, 1);
            BN(0, L) = pn_pxy(0, J) * gradient(1, 0);
            BN(1, L) = pn_pxy(1, J) * gradient(1, 1);
            BN(2, K) = pn_pxy(0, J) * gradient(0, 1) + pn_pxy(1, J) * gradient(0, 0);
            BN(2, L) = pn_pxy(0, J) * gradient(1, 1) + pn_pxy(1, J) * gradient(1, 0);
        }

        if(m_material->update_trial_status(tensor::strain::to_voigt(tensor::strain::to_green(gradient))) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        auto& t_stress = m_material->get_trial_stress();

        const auto sigma = tensor::stress::to_tensor(t_stress);
        for(unsigned J = 0, L = 0, N = 1; J < m_node; ++J, L += m_dof, N += m_dof) {
            const vec t_vec = sigma * pn_pxy.col(J);
            trial_geometry(N, N) = trial_geometry(L, L) = area * thickness * dot(pn_pxy.col(J), t_vec);
            for(auto K = J + 1, M = L + m_dof, P = M + 1; K < m_node; ++K, M += m_dof, P += m_dof) trial_geometry(M, L) = trial_geometry(P, N) = trial_geometry(L, M) = trial_geometry(N, P) = area * thickness * dot(pn_pxy.col(K), t_vec);
        }

        trial_stiffness = area * thickness * BN.t() * m_material->get_trial_stiffness() * BN;
        trial_resistance = area * thickness * BN.t() * t_stress;
    }
    else {
        vec t_strain(3, fill::zeros);
        for(unsigned I = 0; I < m_node; ++I) {
            const auto& t_disp = node_ptr[I].lock()->get_trial_displacement();
            t_strain(0) += t_disp(0) * pn_pxy(0, I);
            t_strain(1) += t_disp(1) * pn_pxy(1, I);
            t_strain(2) += t_disp(0) * pn_pxy(1, I) + t_disp(1) * pn_pxy(0, I);
        }

        if(m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        // trial_stiffness = area * thickness * strain_mat.t() * m_material->get_trial_stiffness() * strain_mat;
        stack_stiffness(trial_stiffness, m_material->get_trial_stiffness(), strain_mat, area * thickness);
        trial_resistance = area * thickness * strain_mat.t() * m_material->get_trial_stress();
    }

    return SUANPAN_SUCCESS;
}

int CP3::commit_status() { return m_material->commit_status(); }

int CP3::clear_status() { return m_material->clear_status(); }

int CP3::reset_status() { return m_material->reset_status(); }

vector<vec> CP3::record(const OutputType T) { return m_material->record(T); }

void CP3::print() {
    suanpan_info("CP3 element connects nodes:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    m_material->print();
}

#ifdef SUANPAN_VTK
#include <vtkTriangle.h>

void CP3::Setup() {
    vtk_cell = vtkSmartPointer<vtkTriangle>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void CP3::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration(), m_dof, m_node);
    else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity(), m_dof, m_node);
    else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

mat CP3::GetData(const OutputType P) {
    vec t_stress(6, fill::zeros);
    if(const auto t_data = m_material->record(P); !t_data.empty()) t_stress(uvec{0, 1, 3}) = t_data[0];
    return repmat(t_stress, 1, m_node);
}

void CP3::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
