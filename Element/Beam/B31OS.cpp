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

#include "B31OS.h"
#include <Domain/DomainBase.h>
#include <Section/Section.h>
#include <Toolbox/IntegrationPlan.h>

B31OS::IntegrationPoint::IntegrationPoint(const double C, const double W, const double L, unique_ptr<Section>&& M)
    : coor(C)
    , weight(W)
    , b_section(std::move(M)) {
    const auto xp = coor + 1., xm = coor - 1.;
    const auto x3p = 3. * coor + 1., x3m = 3. * coor - 1.;

    // nv1=nw1=nf2
    // nv2=nw2=nf4

    const auto dnu = 1. / L;

    const auto nf1 = .25 * (2. + coor) * xm * xm;
    const auto nf2 = .125 * L * xm * xm * xp;
    const auto nf3 = .25 * (2. - coor) * xp * xp;
    const auto nf4 = .125 * L * xm * xp * xp;

    const auto dnf1 = 1.5 * xm * xp / L;
    const auto dnf2 = .25 * xm * x3p;
    const auto dnf3 = -dnf1;
    const auto dnf4 = .25 * xp * x3m;

    const auto ddnf1 = 6. * coor / L / L;
    const auto ddnf2 = x3m / L;
    const auto ddnf3 = -ddnf1;
    const auto ddnf4 = x3p / L;

    strain_mat.zeros(static_cast<unsigned>(SectionType::OS3D), 9);
    // u'
    strain_mat(0, 0) = dnu;
    // v'
    strain_mat(1, 1) = dnf2;
    strain_mat(1, 2) = dnf4;
    // w'
    strain_mat(2, 3) = dnf2;
    strain_mat(2, 4) = dnf4;
    // v''
    strain_mat(3, 1) = ddnf2;
    strain_mat(3, 2) = ddnf4;
    // w''
    strain_mat(4, 3) = ddnf2;
    strain_mat(4, 4) = ddnf4;
    // f
    strain_mat(5, 5) = nf1;
    strain_mat(5, 7) = nf2;
    strain_mat(5, 6) = nf3;
    strain_mat(5, 8) = nf4;
    // f'
    strain_mat(6, 5) = dnf1;
    strain_mat(6, 7) = dnf2;
    strain_mat(6, 6) = dnf3;
    strain_mat(6, 8) = dnf4;
    // f''
    strain_mat(7, 5) = ddnf1;
    strain_mat(7, 7) = ddnf2;
    strain_mat(7, 6) = ddnf3;
    strain_mat(7, 8) = ddnf4;
    // theta_zi, theta_zj, theta_yi, theta_yj
    strain_mat(8, 1) = strain_mat(9, 2) = strain_mat(10, 3) = strain_mat(11, 4) = 1.;
}

B31OS::B31OS(const unsigned T, uvec&& N, const unsigned S, const unsigned O, const unsigned P, const bool F)
    : SectionOSElement3D(T, b_node, b_dof, std::move(N), uvec{S}, F)
    , orientation_tag(O)
    , int_pt_num(P > 20 ? 20 : P) {}

int B31OS::initialize(const shared_ptr<DomainBase>& D) {
    auto& section_proto = D->get<Section>(section_tag(0));

    if(!D->find_orientation(orientation_tag)) {
        suanpan_warning("Element {} cannot find the assigned transformation {}.\n", get_tag(), orientation_tag);
        return SUANPAN_FAIL;
    }

    b_trans = D->get_orientation(orientation_tag)->get_copy();

    if(b_trans->is_nlgeom() != is_nlgeom()) {
        suanpan_warning("Element {} is assigned with an inconsistent transformation {}.\n", get_tag(), orientation_tag);
        return SUANPAN_FAIL;
    }
    if(OrientationType::B3DOS != b_trans->get_orientation_type()) {
        suanpan_warning("Element {} is assigned with an inconsistent transformation {}, use B3DOSL or B3DOSC only.\n", get_tag(), orientation_tag);
        return SUANPAN_FAIL;
    }

    b_trans->set_element_ptr(this);

    access::rw(length) = b_trans->get_length();

    const IntegrationPlan plan(1, int_pt_num, IntegrationType::LOBATTO);

    const mat& section_stiffness = section_proto->get_initial_stiffness();

    mat local_stiffness(9, 9, fill::zeros);
    int_pt.clear();
    int_pt.reserve(int_pt_num);
    for(unsigned I = 0; I < int_pt_num; ++I) {
        int_pt.emplace_back(plan(I, 0), .5 * plan(I, 1), length, section_proto->get_copy());
        int_pt[I].b_section->set_characteristic_length(int_pt[I].weight * length);
        local_stiffness += int_pt[I].strain_mat.t() * section_stiffness * int_pt[I].strain_mat * int_pt[I].weight * length;
    }

    trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);

    if(const auto linear_density = section_proto->get_linear_density(); linear_density > 0.) trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(linear_density);

    return SUANPAN_SUCCESS;
}

int B31OS::update_status() {
    b_trans->update_status();

    // [0]: uniform axial
    // [1]: strong axis bending near node
    // [2]: strong axis bending far node
    // [3]: weak axis bending near node
    // [4]: weak axis bending far node
    // [5]: torsion near node
    // [6]: torsion far node
    // [7]: warping near node
    // [8]: warping far node
    const auto local_deformation = b_trans->to_local_vec(get_trial_displacement());

    mat local_stiffness(9, 9, fill::zeros), local_geometry(9, 9, fill::zeros);
    vec local_resistance(9, fill::zeros);

    for(const auto& I : int_pt) {
        if(I.b_section->update_trial_status(I.strain_mat * local_deformation) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        local_stiffness += I.strain_mat.t() * I.b_section->get_trial_stiffness() * I.strain_mat * I.weight * length;
        local_geometry += I.strain_mat.t() * I.b_section->get_trial_geometry() * I.strain_mat * I.weight * length;
        local_resistance += I.strain_mat.t() * I.b_section->get_trial_resistance() * I.weight * length;
    }

    trial_resistance = b_trans->to_global_vec(local_resistance);

    if(nlgeom) {
        trial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);
        // the geometry matrix from the section is added to elemental geometry matrix if it is a nonlinear geometry analysis
        trial_geometry = b_trans->to_global_geometry_mat(local_resistance) + b_trans->to_global_stiffness_mat(local_geometry);
    }
    else trial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness + local_geometry);

    return SUANPAN_SUCCESS;
}

int B31OS::commit_status() {
    b_trans->commit_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->commit_status();
    return code;
}

int B31OS::clear_status() {
    b_trans->clear_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->clear_status();
    return code;
}

int B31OS::reset_status() {
    b_trans->reset_status();

    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->reset_status();
    return code;
}

std::vector<vec> B31OS::record(const OutputType P) {
    std::vector<vec> data;
    for(const auto& I : int_pt) append_to(data, I.b_section->record(P));
    return data;
}

void B31OS::print() {
    suanpan_info("A spatial beam element.\n");
    for(const auto& I : int_pt) I.b_section->print();
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

void B31OS::Setup() {
    vtk_cell = vtkSmartPointer<vtkLine>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < b_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

void B31OS::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, b_node, fill::zeros);

    if(OutputType::A == type) t_disp = reshape(get_current_acceleration(), b_dof, b_node);
    else if(OutputType::V == type) t_disp = reshape(get_current_velocity(), b_dof, b_node);
    else if(OutputType::U == type) t_disp = reshape(get_current_displacement(), b_dof, b_node);

    t_disp = t_disp.head_rows(6);

    for(unsigned I = 0; I < b_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void B31OS::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(3) + amplifier * mat(reshape(get_current_displacement(), b_dof, b_node)).rows(0, 2).t();
    for(unsigned I = 0; I < b_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
