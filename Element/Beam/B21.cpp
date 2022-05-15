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

#include "B21.h"
#include <Domain/DomainBase.h>
#include <Element/Utility/B2DC.h>
#include <Section/Section.h>
#include <Toolbox/IntegrationPlan.h>

B21::IntegrationPoint::IntegrationPoint(const double C, const double W, unique_ptr<Section>&& M)
    : coor(C)
    , weight(W)
    , b_section(std::forward<unique_ptr<Section>>(M))
    , strain_mat(2, 3, fill::zeros) {}

B21::B21(const unsigned T, uvec&& N, const unsigned S, const unsigned P, const bool F)
    : SectionElement2D(T, b_node, b_dof, std::forward<uvec>(N), uvec{S}, F)
    , int_pt_num(P > 20 ? 20 : P)
    , b_trans(F ? make_unique<B2DC>() : make_unique<B2DL>()) {}

int B21::initialize(const shared_ptr<DomainBase>& D) {
    auto& section_proto = D->get<Section>(section_tag(0));

    b_trans->set_element_ptr(this);

    access::rw(length) = b_trans->get_length();

    const IntegrationPlan plan(1, int_pt_num, IntegrationType::LOBATTO);

    mat local_stiffness(3, 3, fill::zeros);
    int_pt.clear();
    int_pt.reserve(int_pt_num);
    for(unsigned I = 0; I < int_pt_num; ++I) {
        int_pt.emplace_back(plan(I, 0), .5 * plan(I, 1), section_proto->get_copy());
        int_pt[I].strain_mat(0, 0) = 1.;
        int_pt[I].strain_mat(1, 1) = 3. * plan(I, 0) - 1.;
        int_pt[I].strain_mat(1, 2) = 3. * plan(I, 0) + 1.;
        local_stiffness += int_pt[I].strain_mat.t() * int_pt[I].b_section->get_initial_stiffness() * int_pt[I].strain_mat * int_pt[I].weight / length;
    }

    trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);

    if(const auto linear_density = section_proto->get_parameter(ParameterType::LINEARDENSITY); linear_density > 0.) trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(linear_density);

    return SUANPAN_SUCCESS;
}

int B21::update_status() {
    b_trans->update_status();

    // const auto length = b_trans->get_length();

    const auto local_deformation = b_trans->to_local_vec(get_trial_displacement());

    mat local_stiffness(3, 3, fill::zeros);
    vec local_resistance(3, fill::zeros);
    for(const auto& I : int_pt) {
        if(SUANPAN_SUCCESS != I.b_section->update_trial_status(I.strain_mat * local_deformation / length)) return SUANPAN_FAIL;
        local_stiffness += I.strain_mat.t() * I.b_section->get_trial_stiffness() * I.strain_mat * I.weight / length;
        local_resistance += I.strain_mat.t() * I.b_section->get_trial_resistance() * I.weight;
    }

    trial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);
    trial_resistance = b_trans->to_global_vec(local_resistance);

    if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(local_resistance);

    return SUANPAN_SUCCESS;
}

int B21::commit_status() {
    b_trans->commit_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->commit_status();
    return code;
}

int B21::clear_status() {
    b_trans->clear_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->clear_status();
    return code;
}

int B21::reset_status() {
    b_trans->reset_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->reset_status();
    return code;
}

vector<vec> B21::record(const OutputType P) {
    vector<vec> output;
    for(const auto& I : int_pt) for(const auto& J : I.b_section->record(P)) output.emplace_back(J);
    return output;
}

void B21::print() {
    suanpan_info("A classic 2D displacement based beam element using Hermite interpolation functions%s", nlgeom ? " and corotational formulation.\n" : ".\n");
    node_encoding.t().print("The element connects nodes:");
    if(!is_initialized()) return;
    suanpan_info("Section:\n");
    auto J = 1;
    for(const auto& I : int_pt) {
        suanpan_info("IP %d: ", J++);
        I.b_section->print();
    }
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

void B21::Setup() {
    vtk_cell = vtkSmartPointer<vtkLine>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < b_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void B21::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, b_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_acceleration(), b_dof, b_node);
    else if(OutputType::V == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_velocity(), b_dof, b_node);
    else if(OutputType::U == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_displacement(), b_dof, b_node);

    for(unsigned I = 0; I < b_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void B21::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * mat(reshape(get_current_displacement(), b_dof, b_node).t()).cols(0, 1);
    for(unsigned I = 0; I < b_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
