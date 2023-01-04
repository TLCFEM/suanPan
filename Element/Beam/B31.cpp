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

#include "B31.h"
#include <Domain/DomainBase.h>
#include <Section/Section.h>
#include <Toolbox/IntegrationPlan.h>

const span B31::b_span(0, 2);

B31::IntegrationPoint::IntegrationPoint(const double C, const double W, unique_ptr<Section>&& M)
    : coor(C)
    , weight(W)
    , b_section(std::forward<unique_ptr<Section>>(M))
    , strain_mat(3, 6, fill::zeros) {}

B31::B31(const unsigned T, uvec&& N, const unsigned S, const unsigned O, const unsigned P, const bool F)
    : SectionElement3D(T, b_node, b_dof, std::forward<uvec>(N), uvec{S}, F)
    , orientation_tag(O)
    , int_pt_num(P > 20 ? 20 : P) {}

int B31::initialize(const shared_ptr<DomainBase>& D) {
    auto& sec_proto = D->get<Section>(section_tag(0));

    const mat sec_stiff = sec_proto->get_initial_stiffness()(b_span, b_span);

    if(!D->find_orientation(orientation_tag)) {
        suanpan_warning("Element {} cannot find the assigned transformation {}.\n", get_tag(), orientation_tag);
        return SUANPAN_FAIL;
    }

    b_trans = D->get_orientation(orientation_tag)->get_copy();

    if(b_trans->is_nlgeom() != is_nlgeom()) {
        suanpan_warning("Element {} is assigned with an inconsistent transformation {}.\n", get_tag(), orientation_tag);
        return SUANPAN_FAIL;
    }

    b_trans->set_element_ptr(this);

    access::rw(length) = b_trans->get_length();

    const IntegrationPlan plan(1, int_pt_num, IntegrationType::LOBATTO);

    mat local_stiffness(6, 6, fill::zeros);
    int_pt.clear();
    int_pt.reserve(int_pt_num);
    for(unsigned I = 0; I < int_pt_num; ++I) {
        int_pt.emplace_back(plan(I, 0), .5 * plan(I, 1), sec_proto->get_copy());
        int_pt[I].strain_mat(0, 0) = 1.;
        int_pt[I].strain_mat(1, 1) = int_pt[I].strain_mat(2, 3) = 3. * plan(I, 0) - 1.;
        int_pt[I].strain_mat(1, 2) = int_pt[I].strain_mat(2, 4) = 3. * plan(I, 0) + 1.;
        local_stiffness += int_pt[I].strain_mat.t() * sec_stiff * int_pt[I].strain_mat * int_pt[I].weight / length;
    }
    access::rw(torsion_stiff) = 1E3 * local_stiffness.max();
    local_stiffness(5, 5) = torsion_stiff;

    trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);

    if(const auto linear_density = sec_proto->get_parameter(ParameterType::LINEARDENSITY); linear_density > 0.) trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(linear_density);

    return SUANPAN_SUCCESS;
}

int B31::update_status() {
    b_trans->update_status();

    const auto local_deformation = b_trans->to_local_vec(get_trial_displacement());

    mat local_stiffness(6, 6, fill::zeros);
    vec local_resistance(6, fill::zeros);
    for(const auto& I : int_pt) {
        if(I.b_section->update_trial_status(I.strain_mat * local_deformation / length) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        local_stiffness += I.strain_mat.t() * I.b_section->get_trial_stiffness()(b_span, b_span) * I.strain_mat * I.weight / length;
        local_resistance += I.strain_mat.t() * I.b_section->get_trial_resistance()(b_span) * I.weight;
    }
    local_resistance(5) = (local_stiffness(5, 5) = torsion_stiff) * local_deformation(5);

    trial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);
    trial_resistance = b_trans->to_global_vec(local_resistance);

    if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(local_resistance);

    return SUANPAN_SUCCESS;
}

int B31::commit_status() {
    b_trans->commit_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->commit_status();
    return code;
}

int B31::clear_status() {
    b_trans->clear_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->clear_status();
    return code;
}

int B31::reset_status() {
    b_trans->reset_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->reset_status();
    return code;
}

vector<vec> B31::record(const OutputType P) {
    vector<vec> output;
    for(const auto& I : int_pt) for(const auto& J : I.b_section->record(P)) output.emplace_back(J);
    return output;
}

void B31::print() {
    suanpan_info("A spatial beam element.\n");
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

void B31::Setup() {
    vtk_cell = vtkSmartPointer<vtkLine>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < b_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

void B31::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, b_node, fill::zeros);

    if(OutputType::A == type) t_disp = reshape(get_current_acceleration(), b_dof, b_node);
    else if(OutputType::V == type) t_disp = reshape(get_current_velocity(), b_dof, b_node);
    else if(OutputType::U == type) t_disp = reshape(get_current_displacement(), b_dof, b_node);

    for(unsigned I = 0; I < b_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void B31::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(3) + amplifier * mat(reshape(get_current_displacement(), b_dof, b_node)).rows(0, 2).t();
    for(unsigned I = 0; I < b_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
