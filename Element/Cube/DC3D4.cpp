/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "DC3D4.h"
#include <Domain/DOF.h>
#include <Domain/DomainBase.h>
#include <Material/Material3D/Material3D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/shapeFunction.h>

const uvec DC3D4::u_dof{0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14};
const uvec DC3D4::d_dof{3, 7, 11, 15};

DC3D4::DC3D4(const unsigned T, uvec&& N, const unsigned M, const double CL, const double RR)
    : MaterialElement3D(T, c_node, c_dof, std::forward<uvec>(N), uvec{M}, false, {DOF::X, DOF::Y, DOF::Z, DOF::DMG})
    , release_rate(RR) { access::rw(characteristic_length) = CL; }

int DC3D4::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    c_material = material_proto->get_copy();

    mat ele_coor(c_node, c_node);
    ele_coor.col(0).fill(1.);
    ele_coor.cols(1, 3) = get_coordinate(3);

    access::rw(volume) = det(ele_coor) / 6.;

    if(0. >= characteristic_length) access::rw(characteristic_length) = 2. * pow(volume, 1. / 3.);

    const mat inv_coor = inv(ele_coor);

    n_mat = mean(ele_coor) * inv_coor;

    pn_mat = inv_coor.rows(1, 3);

    b_mat.zeros(6, 12);
    for(unsigned J = 0, K = 0, L = 1, M = 2; J < c_node; ++J, K += 3, L += 3, M += 3) {
        b_mat(0, K) = b_mat(3, L) = b_mat(5, M) = pn_mat(0, J);
        b_mat(3, K) = b_mat(1, L) = b_mat(4, M) = pn_mat(1, J);
        b_mat(5, K) = b_mat(4, L) = b_mat(2, M) = pn_mat(2, J);
    }
    initial_stiffness.zeros(c_size, c_size);
    initial_stiffness(u_dof, u_dof) = b_mat.t() * c_material->get_initial_stiffness() * b_mat;
    initial_stiffness(d_dof, d_dof) = release_rate / characteristic_length * n_mat.t() * n_mat + release_rate * characteristic_length * pn_mat.t() * pn_mat;
    trial_stiffness = current_stiffness = initial_stiffness *= volume;

    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

int DC3D4::update_status() {
    const auto t_disp = get_trial_displacement();
    const vec t_damage = t_disp(d_dof);

    if(c_material->update_trial_status(b_mat * t_disp(u_dof)) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    const auto pow_term = 1. - dot(t_damage, n_mat);
    const auto damage = pow(pow_term, 2.);

    trial_stiffness.zeros(c_size, c_size);
    trial_resistance.zeros(c_size);

    trial_stiffness(u_dof, u_dof) = damage * b_mat.t() * c_material->get_trial_stiffness() * b_mat;
    trial_stiffness(u_dof, d_dof) = -2. * pow_term * b_mat.t() * c_material->get_trial_stress() * n_mat;
    trial_stiffness(d_dof, d_dof) = n_mat.t() * n_mat * (2. * maximum_energy + release_rate / characteristic_length) + release_rate * characteristic_length * pn_mat.t() * pn_mat;

    trial_resistance(u_dof) = damage * b_mat.t() * c_material->get_trial_stress();
    trial_resistance(d_dof) = trial_stiffness(d_dof, d_dof) * t_damage;
    trial_resistance(d_dof) -= 2. * n_mat.t() * maximum_energy;

    trial_stiffness *= volume;
    trial_resistance *= volume;

    return SUANPAN_SUCCESS;
}

int DC3D4::commit_status() {
    c_phase.commit_status(c_material);
    maximum_energy = std::max(maximum_energy, c_phase.strain_energy);

    return c_material->commit_status();
}

int DC3D4::clear_status() {
    c_phase.clear_status();
    maximum_energy = 0.;

    return c_material->clear_status();
}

int DC3D4::reset_status() { return c_material->reset_status(); }

vector<vec> DC3D4::record(const OutputType T) {
    if(T == OutputType::DAMAGE) return {get_current_displacement()(d_dof)};

    return c_material->record(T);
}

void DC3D4::print() {
    node_encoding.t().print("DC3D4 element connects:");
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    c_material->print();
    suanpan_info("Strain:\t");
    c_material->get_current_strain().t().print();
    suanpan_info("Stress:\t");
    c_material->get_current_stress().t().print();
}

#ifdef SUANPAN_VTK
#include <vtkTetra.h>

void DC3D4::Setup() {
    vtk_cell = vtkSmartPointer<vtkTetra>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < c_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

mat DC3D4::GetData(const OutputType P) {
    if(P == OutputType::DAMAGE) {
        mat t_damage(6, c_node, fill::zeros);
        t_damage.row(0) = get_current_displacement()(d_dof).t();
        return t_damage;
    }

    return {};
}

void DC3D4::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, c_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 2) = reshape(get_current_acceleration()(u_dof), 3, c_node);
    else if(OutputType::V == type) t_disp.rows(0, 2) = reshape(get_current_velocity()(u_dof), 3, c_node);
    else if(OutputType::U == type) t_disp.rows(0, 2) = reshape(get_current_displacement()(u_dof), 3, c_node);

    for(unsigned I = 0; I < c_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void DC3D4::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(3) + amplifier * reshape(get_current_displacement()(u_dof), 3, c_node).t();
    for(unsigned I = 0; I < c_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
