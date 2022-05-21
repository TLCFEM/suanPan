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

#include "F31.h"
#include <Domain/DomainBase.h>
#include <Section/Section.h>
#include <Toolbox/IntegrationPlan.h>

const span F31::b_span(0, 2);

F31::IntegrationPoint::IntegrationPoint(const double C, const double W, unique_ptr<Section>&& M)
    : coor(C)
    , weight(W)
    , b_section(std::forward<unique_ptr<Section>>(M))
    , strain_mat(3, 6, fill::zeros) {}

F31::F31(const unsigned T, uvec&& N, const unsigned S, const unsigned O, const unsigned P, const bool F)
    : SectionElement3D(T, b_node, b_dof, std::forward<uvec>(N), uvec{S}, F)
    , int_pt_num(P > 20 ? 20 : P)
    , orientation_tag(O) {}

int F31::initialize(const shared_ptr<DomainBase>& D) {
    auto& sec_proto = D->get<Section>(section_tag(0));

    if(!D->find_orientation(orientation_tag)) {
        suanpan_warning("Element %u cannot find the assigned transformation.\n", get_tag());
        return SUANPAN_FAIL;
    }

    b_trans = D->get_orientation(orientation_tag)->get_copy();

    if(b_trans->is_nlgeom() != is_nlgeom()) {
        suanpan_warning("Element %u is assigned with an inconsistent transformation.\n", get_tag());
        return SUANPAN_FAIL;
    }

    b_trans->set_element_ptr(this);

    access::rw(length) = b_trans->get_length();

    const mat sec_stiff = sec_proto->get_initial_stiffness()(b_span, b_span);

    const IntegrationPlan plan(1, int_pt_num, IntegrationType::LOBATTO);

    initial_local_flexibility.zeros(6, 6);
    int_pt.clear();
    int_pt.reserve(int_pt_num);
    for(unsigned I = 0; I < int_pt_num; ++I) {
        int_pt.emplace_back(plan(I, 0), .5 * plan(I, 1), sec_proto->get_copy());
        int_pt[I].strain_mat(0, 0) = 1.;
        int_pt[I].strain_mat(1, 1) = int_pt[I].strain_mat(2, 3) = .5 * plan(I, 0) - .5;
        int_pt[I].strain_mat(1, 2) = int_pt[I].strain_mat(2, 4) = .5 * plan(I, 0) + .5;
        // factor .5 moved to weight
        initial_local_flexibility += int_pt[I].strain_mat.t() * solve(sec_stiff, int_pt[I].strain_mat * int_pt[I].weight * length);
    }
    access::rw(torsion_stiff) = 1E-3 * vec(initial_local_flexibility.diag()).head(5).min();
    initial_local_flexibility(5, 5) = torsion_stiff;
    trial_local_flexibility = current_local_flexibility = initial_local_flexibility;

    trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(inv(initial_local_flexibility));

    if(const auto linear_density = sec_proto->get_parameter(ParameterType::LINEARDENSITY); linear_density > 0.) trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(linear_density);

    trial_local_deformation = current_local_deformation.zeros(6);
    trial_local_resistance = current_local_resistance.zeros(6);

    return SUANPAN_SUCCESS;
}

int F31::update_status() {
    b_trans->update_status();

    vec residual_deformation = -trial_local_deformation;
    // transform global deformation to local one (remove rigid body motion)
    trial_local_deformation = b_trans->to_local_vec(get_trial_displacement());
    // initial residual be aware of how to compute it
    residual_deformation += trial_local_deformation;

    trial_local_resistance += solve(trial_local_flexibility, residual_deformation);

    residual_deformation.zeros();
    trial_local_flexibility.zeros();
    for(const auto& I : int_pt) {
        const vec target_resistance = I.strain_mat * trial_local_resistance;
        // update status
        vec incre_deformation;
        if(!solve(incre_deformation, I.b_section->get_trial_stiffness()(b_span, b_span), target_resistance - I.b_section->get_trial_resistance()(b_span))) return SUANPAN_FAIL;
        if(SUANPAN_SUCCESS != I.b_section->update_trial_status(I.b_section->get_trial_deformation() + incre_deformation)) return SUANPAN_FAIL;
        // collect new flexibility and deformation
        trial_local_flexibility += I.weight * length * I.strain_mat.t() * solve(I.b_section->get_trial_stiffness()(b_span, b_span), I.strain_mat);
        residual_deformation += I.weight * length * I.strain_mat.t() * solve(I.b_section->get_trial_stiffness()(b_span, b_span), target_resistance - I.b_section->get_trial_resistance()(b_span));
    }
    trial_local_flexibility(5, 5) = torsion_stiff;

    trial_local_resistance -= solve(trial_local_flexibility, residual_deformation);

    trial_stiffness = b_trans->to_global_stiffness_mat(inv(trial_local_flexibility));
    trial_resistance = b_trans->to_global_vec(trial_local_resistance);

    if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(trial_local_resistance);

    return SUANPAN_SUCCESS;
}

int F31::clear_status() {
    b_trans->clear_status();

    trial_local_flexibility = current_local_flexibility = initial_local_flexibility;
    trial_local_deformation = current_local_deformation.zeros();
    trial_local_resistance = current_local_resistance.zeros();
    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->clear_status();
    return code;
}

int F31::commit_status() {
    b_trans->commit_status();

    current_local_flexibility = trial_local_flexibility;
    current_local_deformation = trial_local_deformation;
    current_local_resistance = trial_local_resistance;
    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->commit_status();
    return code;
}

int F31::reset_status() {
    b_trans->reset_status();

    trial_local_flexibility = current_local_flexibility;
    trial_local_deformation = current_local_deformation;
    trial_local_resistance = current_local_resistance;
    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->reset_status();
    return code;
}

vector<vec> F31::record(const OutputType P) {
    vector<vec> output;
    for(const auto& I : int_pt) for(const auto& J : I.b_section->record(P)) output.emplace_back(J);
    return output;
}

void F31::print() {
    suanpan_info("A 2D force based beam element%s.\n", nlgeom ? " and corotational formulation" : "");
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

void F31::Setup() {
    vtk_cell = vtkSmartPointer<vtkLine>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < b_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

void F31::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, b_node, fill::zeros);

    if(OutputType::A == type) t_disp = reshape(get_current_acceleration(), b_dof, b_node);
    else if(OutputType::V == type) t_disp = reshape(get_current_velocity(), b_dof, b_node);
    else if(OutputType::U == type) t_disp = reshape(get_current_displacement(), b_dof, b_node);

    for(unsigned I = 0; I < b_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void F31::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(3) + amplifier * mat(reshape(get_current_displacement(), b_dof, b_node)).rows(0, 2).t();
    for(unsigned I = 0; I < b_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
