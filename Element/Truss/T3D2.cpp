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

#include "T3D2.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

T3D2::T3D2(const unsigned T, uvec&& N, const unsigned M, const double A, const bool F, const bool UA, const bool LS)
    : MaterialElement1D(T, t_node, t_dof, std::forward<uvec>(N), uvec{M}, F, {DOF::U1, DOF::U2, DOF::U3})
    , area(A)
    , t_trans(F ? make_unique<T3DC>() : make_unique<T3DL>())
    , update_area(UA)
    , log_strain(LS) {}

int T3D2::initialize(const shared_ptr<DomainBase>& D) {
    t_trans->set_element_ptr(this);

    access::rw(length) = t_trans->get_length();

    t_material = D->get<Material>(material_tag(0))->get_copy();

    trial_stiffness = current_stiffness = initial_stiffness = t_trans->to_global_stiffness_mat(area / length * t_material->get_initial_stiffness());

    if(const auto t_density = t_material->get_density(); t_density > 0.) trial_mass = current_mass = initial_mass = t_trans->to_global_mass_mat(t_density * area);

    return SUANPAN_SUCCESS;
}

int T3D2::update_status() {
    auto new_area = area;

    t_trans->update_status();

    if(nlgeom) {
        const auto new_length = t_trans->get_length();

        double d_strain;

        if(t_material->update_trial_status(log_strain ? log((d_strain = new_length) / length) : new_length / (d_strain = length) - 1.) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        const auto t_stiff = update_area ? -t_material->get_trial_stress().at(0) * (new_area *= length / new_length) / new_length : 0.;

        trial_stiffness = t_trans->to_global_stiffness_mat(t_material->get_trial_stiffness() * new_area / d_strain + t_stiff);
        trial_geometry = t_trans->to_global_geometry_mat(new_area / new_length * t_material->get_trial_stress());
    }
    else {
        if(t_material->update_trial_status(t_trans->to_local_vec(get_trial_displacement()) / length) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_stiffness = t_trans->to_global_stiffness_mat(area / length * t_material->get_trial_stiffness());
    }

    trial_resistance = t_trans->to_global_vec(new_area * t_material->get_trial_stress());

    return SUANPAN_SUCCESS;
}

int T3D2::commit_status() {
    t_trans->commit_status();

    return t_material->commit_status();
}

int T3D2::clear_status() {
    t_trans->clear_status();

    return t_material->clear_status();
}

int T3D2::reset_status() {
    t_trans->reset_status();

    return t_material->reset_status();
}

vector<vec> T3D2::record(const OutputType T) { return t_material->record(T); }

void T3D2::print() {
    suanpan_info("A 3D truss element with ");
    if(nlgeom)
        suanpan_info("corotational formulation, assuming constant {} and {} strain.", update_area ? "volume" : "area", log_strain ? "logarithmic" : "engineering");
    else
        suanpan_info("linear formulation.");
    suanpan_info(" The nodes connected are:", node_encoding);
    suanpan_info("The area is {:.4E}. The initial element length is {:.4E}.\n", area, length);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    t_material->print();
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

void T3D2::Setup() {
    vtk_cell = vtkSmartPointer<vtkLine>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < t_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

void T3D2::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, t_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 2) = reshape(get_current_acceleration(), t_dof, t_node);
    else if(OutputType::V == type) t_disp.rows(0, 2) = reshape(get_current_velocity(), t_dof, t_node);
    else if(OutputType::U == type) t_disp.rows(0, 2) = reshape(get_current_displacement(), t_dof, t_node);

    for(unsigned I = 0; I < t_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void T3D2::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(3) + amplifier * reshape(get_current_displacement(), t_dof, t_node).t();
    for(unsigned I = 0; I < t_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
