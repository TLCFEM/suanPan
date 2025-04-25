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

#include "B21H.h"

#include <Domain/DomainBase.h>
#include <Section/Section.h>
#include <Toolbox/IntegrationPlan.h>

B21H::IntegrationPoint::IntegrationPoint(const double C, const double W, unique_ptr<Section>&& M)
    : coor(C)
    , weight(W)
    , b_section(std::move(M))
    , strain_mat(2, 3, fill::zeros) {
    strain_mat(0, 0) = 1.;
    strain_mat(1, 1) = 3. * coor - 1.;
    strain_mat(1, 2) = 3. * coor + 1.;
}

B21H::B21H(const unsigned T, uvec&& N, const unsigned S, const double L, const bool F)
    : SectionElement2D(T, b_node, b_dof, std::move(N), uvec{S}, F)
    , hinge_length(L > .5 ? .5 : L)
    , b_trans(F ? make_unique<B2DC>() : make_unique<B2DL>()) {}

int B21H::initialize(const shared_ptr<DomainBase>& D) {
    auto& section_proto = D->get<Section>(section_tag(0));

    b_trans->set_element_ptr(this);

    access::rw(length) = b_trans->get_length();

    const auto elastic_length = 1. - 2. * hinge_length;

    // build up the elastic interior
    const IntegrationPlan elastic_plan(1, 2, IntegrationType::GAUSS);
    elastic_int_pt.clear();
    elastic_int_pt.reserve(elastic_plan.n_rows);
    for(unsigned I = 0; I < elastic_plan.n_rows; ++I) elastic_int_pt.emplace_back(elastic_plan(I, 0) * elastic_length, elastic_plan(I, 1) * elastic_length / 2., section_proto->get_copy());

    int_pt.clear();
    int_pt.reserve(4);
    int_pt.emplace_back(-1., .25 * hinge_length, section_proto->get_copy());
    int_pt.emplace_back(4. / 3. * hinge_length - 1., .75 * hinge_length, section_proto->get_copy());
    int_pt.emplace_back(1. - 4. / 3. * hinge_length, .75 * hinge_length, section_proto->get_copy());
    int_pt.emplace_back(1., .25 * hinge_length, section_proto->get_copy());

    const auto& elastic_section_stiffness = section_proto->get_initial_stiffness();
    // elastic part will be reused in computation
    elastic_local_stiffness.zeros(3, 3);
    for(auto& I : elastic_int_pt) elastic_local_stiffness += I.strain_mat.t() * elastic_section_stiffness * I.strain_mat * I.weight / length;

    auto local_stiffness = elastic_local_stiffness;
    for(auto& I : int_pt) local_stiffness += I.strain_mat.t() * I.b_section->get_initial_stiffness() * I.strain_mat * I.weight / length;
    trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);

    if(const auto linear_density = section_proto->get_linear_density(); linear_density > 0.) trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(linear_density);

    return SUANPAN_SUCCESS;
}

int B21H::update_status() {
    b_trans->update_status();

    const auto local_deformation = b_trans->to_local_vec(get_trial_displacement());

    auto local_stiffness = elastic_local_stiffness;
    vec local_resistance = elastic_local_stiffness * local_deformation;
    for(const auto& I : int_pt) {
        if(SUANPAN_SUCCESS != I.b_section->update_trial_status(I.strain_mat * local_deformation / length)) return SUANPAN_FAIL;
        local_stiffness += I.strain_mat.t() * I.b_section->get_trial_stiffness() * I.strain_mat * I.weight / length;
        local_resistance += I.strain_mat.t() * I.b_section->get_trial_resistance() * I.weight;
    }

    for(const auto& I : elastic_int_pt) I.b_section->update_trial_status(I.strain_mat * local_deformation / length);

    trial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);
    trial_resistance = b_trans->to_global_vec(local_resistance);

    if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(local_resistance);

    return SUANPAN_SUCCESS;
}

int B21H::commit_status() {
    b_trans->commit_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->commit_status();
    return code;
}

int B21H::clear_status() {
    b_trans->clear_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->clear_status();
    return code;
}

int B21H::reset_status() {
    b_trans->reset_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->reset_status();
    return code;
}

std::vector<vec> B21H::record(const OutputType P) {
    std::vector<vec> data;
    append_to(data, int_pt[0].b_section->record(P));
    append_to(data, int_pt[1].b_section->record(P));
    for(const auto& I : elastic_int_pt) append_to(data, I.b_section->record(P));
    append_to(data, int_pt[2].b_section->record(P));
    append_to(data, int_pt[3].b_section->record(P));
    return data;
}

void B21H::print() {
    suanpan_info("A 2D beam element with lumped end plasticity (hinges){}", nlgeom ? " and corotational formulation.\n" : ".\n");
    suanpan_info("The plastic hinge length is: {:.3f}.\n", hinge_length);
    suanpan_info("The element connects nodes:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Section:\n");
    auto J = 1;
    for(const auto& I : int_pt) {
        suanpan_info("IP {}: ", J++);
        I.b_section->print();
    }
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

void B21H::Setup() {
    vtk_cell = vtkSmartPointer<vtkLine>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < b_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void B21H::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, b_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_acceleration(), b_dof, b_node);
    else if(OutputType::V == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_velocity(), b_dof, b_node);
    else if(OutputType::U == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_displacement(), b_dof, b_node);

    for(unsigned I = 0; I < b_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void B21H::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * mat(reshape(get_current_displacement(), b_dof, b_node).t()).cols(0, 1);
    for(unsigned I = 0; I < b_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
