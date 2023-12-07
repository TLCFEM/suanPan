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

#include "DCP3.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/utility.h>

const uvec DCP3::u_dof{0, 1, 3, 4, 6, 7};
const uvec DCP3::d_dof{2, 5, 8};

DCP3::DCP3(const unsigned T, uvec&& NT, const unsigned MT, const double CL, const double RR, const double TH)
    : MaterialElement2D(T, m_node, m_dof, std::move(NT), uvec{MT}, false, {DOF::U1, DOF::U2, DOF::DMG})
    , release_rate(RR)
    , thickness(TH) { access::rw(characteristic_length) = CL; }

int DCP3::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(PlaneType::E == material_proto->get_plane_type()) suanpan::hacker(thickness) = 1.;

    m_material = material_proto->get_copy();

    mat ele_coor(m_node, m_node);
    ele_coor.col(0).fill(1.);
    ele_coor.cols(1, 2) = get_coordinate(2);

    access::rw(area) = .5 * det(ele_coor);

    if(0. >= characteristic_length) access::rw(characteristic_length) = 2. * sqrt(area);

    const mat inv_coor = inv(ele_coor);
    pn_mat = inv_coor.rows(1, 2);

    n_mat = mean(ele_coor) * inv_coor;

    b_mat.zeros(3, 6);
    for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += 2, L += 2) {
        b_mat(0, K) = b_mat(2, L) = pn_mat(0, J);
        b_mat(2, K) = b_mat(1, L) = pn_mat(1, J);
    }
    initial_stiffness.zeros(m_size, m_size);
    initial_stiffness(u_dof, u_dof) = b_mat.t() * m_material->get_initial_stiffness() * b_mat;
    initial_stiffness(d_dof, d_dof) = release_rate / characteristic_length * n_mat.t() * n_mat + release_rate * characteristic_length * pn_mat.t() * pn_mat;
    trial_stiffness = current_stiffness = initial_stiffness *= area * thickness;

    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

int DCP3::update_status() {
    vec t_strain(3, fill::zeros), t_damage(m_node);
    for(unsigned I = 0; I < m_node; ++I) {
        const auto& t_disp = node_ptr[I].lock()->get_trial_displacement();
        t_strain(0) += t_disp(0) * pn_mat(0, I);
        t_strain(1) += t_disp(1) * pn_mat(1, I);
        t_strain(2) += t_disp(0) * pn_mat(1, I) + t_disp(1) * pn_mat(0, I);
        t_damage(I) = t_disp(2);
    }

    if(m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    const auto pow_term = 1. - dot(t_damage, n_mat);
    const auto damage = pow(pow_term, 2.);

    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    trial_stiffness(u_dof, u_dof) = damage * b_mat.t() * m_material->get_trial_stiffness() * b_mat;
    trial_stiffness(u_dof, d_dof) = -2. * pow_term * b_mat.t() * m_material->get_trial_stress() * n_mat;
    trial_stiffness(d_dof, d_dof) = n_mat.t() * n_mat * (2. * maximum_energy + release_rate / characteristic_length) + release_rate * characteristic_length * pn_mat.t() * pn_mat;

    trial_resistance(u_dof) = damage * b_mat.t() * m_material->get_trial_stress();
    trial_resistance(d_dof) = trial_stiffness(d_dof, d_dof) * t_damage;
    trial_resistance(d_dof) -= 2. * n_mat.t() * maximum_energy;

    trial_stiffness *= area * thickness;
    trial_resistance *= area * thickness;

    return SUANPAN_SUCCESS;
}

int DCP3::commit_status() {
    m_phase.commit_status(m_material);
    maximum_energy = std::max(maximum_energy, m_phase.strain_energy);

    return m_material->commit_status();
}

int DCP3::clear_status() {
    m_phase.clear_status();
    maximum_energy = 0.;

    return m_material->clear_status();
}

int DCP3::reset_status() { return m_material->reset_status(); }

vector<vec> DCP3::record(const OutputType P) {
    if(P == OutputType::DAMAGE) return {get_current_displacement()(d_dof)};

    return m_material->record(P);
}

void DCP3::print() {
    suanpan_info("DCP3 element connects nodes:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    m_material->print();
}

#ifdef SUANPAN_VTK
#include <vtkTriangle.h>

void DCP3::Setup() {
    vtk_cell = vtkSmartPointer<vtkTriangle>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void DCP3::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration()(u_dof), 2, m_node);
    else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity()(u_dof), 2, m_node);
    else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement()(u_dof), 2, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

mat DCP3::GetData(const OutputType P) {
    if(P == OutputType::DAMAGE) {
        mat t_damage(6, m_node, fill::zeros);
        for(unsigned I = 0; I < m_node; ++I) t_damage(0, I) = node_ptr[I].lock()->get_current_displacement()(2);
        return t_damage;
    }

    vec t_stress(6, fill::zeros);
    if(const auto t_data = m_material->record(P); !t_data.empty()) t_stress(uvec{0, 1, 3}) = t_data[0];
    return repmat(t_stress, 1, m_node);
}

void DCP3::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement()(u_dof), 2, m_node).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
