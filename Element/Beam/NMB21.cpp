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

#include "NMB21.h"

#include <Domain/DomainBase.h>
#include <Element/Utility/B2DC.h>
#include <Recorder/OutputType.h>
#include <Section/Section.h>

NMB21::NMB21(const unsigned T, uvec&& N, const unsigned S, const bool R)
    : SectionNMElement2D(T, b_node, b_dof, std::move(N), uvec{S}, R)
    , b_trans(R ? make_unique<B2DC>() : make_unique<B2DL>()) {}

int NMB21::initialize(const shared_ptr<DomainBase>& D) {
    b_trans->set_element_ptr(this);

    access::rw(length) = b_trans->get_length();

    b_section = D->get<Section>(section_tag(0))->get_copy();

    trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(b_section->get_initial_stiffness() / length);

    if(const auto linear_density = b_section->get_linear_density(); linear_density > 0.) trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(linear_density);

    return SUANPAN_SUCCESS;
}

int NMB21::update_status() {
    b_trans->update_status();

    if(b_section->update_trial_status(b_trans->to_local_vec(get_trial_displacement()) / length) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    trial_stiffness = b_trans->to_global_stiffness_mat(b_section->get_trial_stiffness() / length);
    trial_resistance = b_trans->to_global_vec(b_section->get_trial_resistance());

    if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(b_section->get_trial_resistance());

    return SUANPAN_SUCCESS;
}

int NMB21::commit_status() {
    b_trans->commit_status();
    return b_section->commit_status();
}

int NMB21::clear_status() {
    b_trans->clear_status();
    return b_section->clear_status();
}

int NMB21::reset_status() {
    b_trans->reset_status();
    return b_section->reset_status();
}

std::vector<vec> NMB21::record(const OutputType P) {
    if(P == OutputType::BEAME) return {b_section->get_current_deformation() * length};
    if(P == OutputType::BEAMS) return {b_section->get_current_resistance()};

    return b_section->record(P);
}

void NMB21::print() {
    suanpan_info("A planar beam element using N-M interaction section.\n");
    if(b_section) b_section->print();
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

void NMB21::Setup() {
    vtk_cell = vtkSmartPointer<vtkLine>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < b_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void NMB21::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, b_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_acceleration(), b_dof, b_node);
    else if(OutputType::V == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_velocity(), b_dof, b_node);
    else if(OutputType::U == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_displacement(), b_dof, b_node);

    for(unsigned I = 0; I < b_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void NMB21::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * mat(reshape(get_current_displacement(), b_dof, b_node)).rows(0, 1).t();
    for(unsigned I = 0; I < b_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
