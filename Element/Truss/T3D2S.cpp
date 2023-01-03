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

#include "T3D2S.h"
#include <Domain/DomainBase.h>
#include <Element/Utility/T3DC.h>
#include <Section/Section.h>

T3D2S::T3D2S(const unsigned T, uvec&& N, const unsigned M, const bool F, const bool LS)
    : SectionElement1D(T, t_node, t_dof, std::forward<uvec>(N), uvec{M}, F, {DOF::U1, DOF::U2, DOF::U3})
    , t_trans(F ? make_unique<T3DC>() : make_unique<T3DL>())
    , log_strain(LS) {}

int T3D2S::initialize(const shared_ptr<DomainBase>& D) {
    t_trans->set_element_ptr(this);

    access::rw(length) = t_trans->get_length();

    t_section = D->get<Section>(section_tag(0))->get_copy();

    trial_stiffness = current_stiffness = initial_stiffness = t_trans->to_global_stiffness_mat(t_section->get_initial_stiffness() / length);

    if(const auto linear_density = t_section->get_parameter(ParameterType::LINEARDENSITY); linear_density > 0.) trial_mass = current_mass = initial_mass = t_trans->to_global_mass_mat(linear_density);

    return SUANPAN_SUCCESS;
}

int T3D2S::update_status() {
    t_trans->update_status();

    if(nlgeom) {
        const auto new_length = t_trans->get_length();

        double d_strain;

        if(t_section->update_trial_status(log_strain ? log((d_strain = new_length) / length) : new_length / (d_strain = length) - 1.) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_stiffness = t_trans->to_global_stiffness_mat(t_section->get_trial_stiffness() / d_strain);
        trial_geometry = t_trans->to_global_geometry_mat(t_section->get_trial_resistance() / new_length);
    }
    else {
        if(t_section->update_trial_status(t_trans->to_local_vec(get_trial_displacement()) / length) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_stiffness = t_trans->to_global_stiffness_mat(t_section->get_trial_stiffness() / length);
    }

    trial_resistance = t_trans->to_global_vec(t_section->get_trial_resistance());

    return SUANPAN_SUCCESS;
}

int T3D2S::commit_status() {
    t_trans->commit_status();

    return t_section->commit_status();
}

int T3D2S::clear_status() {
    t_trans->clear_status();

    return t_section->clear_status();
}

int T3D2S::reset_status() {
    t_trans->reset_status();

    return t_section->reset_status();
}

vector<vec> T3D2S::record(const OutputType T) { return t_section->record(T); }

void T3D2S::print() {
    suanpan_info("A 3D truss element with ");
    if(nlgeom)
        suanpan_info("corotational formulation, assuming constant area and {} strain.", log_strain ? "logarithmic" : "engineering");
    else
        suanpan_info("linear formulation.");
    suanpan_info(" The nodes connected are:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Section:\n");
    t_section->print();
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

void T3D2S::Setup() {
    vtk_cell = vtkSmartPointer<vtkLine>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < t_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

void T3D2S::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, t_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 2) = reshape(get_current_acceleration(), t_dof, t_node);
    else if(OutputType::V == type) t_disp.rows(0, 2) = reshape(get_current_velocity(), t_dof, t_node);
    else if(OutputType::U == type) t_disp.rows(0, 2) = reshape(get_current_displacement(), t_dof, t_node);

    for(unsigned I = 0; I < t_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void T3D2S::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(3) + amplifier * reshape(get_current_displacement(), t_dof, t_node).t();
    for(unsigned I = 0; I < t_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
