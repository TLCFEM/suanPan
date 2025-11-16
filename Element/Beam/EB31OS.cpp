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

#include "EB31OS.h"

#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>

EB31OS::EB31OS(const unsigned T, uvec&& N, vec&& P, const unsigned O, const bool F)
    : SectionOSElement3D(T, b_node, b_dof, std::move(N), uvec{}, F)
    , orientation_tag(O)
    , property(std::move(P)) {}

int EB31OS::initialize(const shared_ptr<DomainBase>& D) {
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

    const auto& E = property(0);
    const auto& G = property(1);
    const auto& A = property(2);
    const auto& IZ = property(3);
    const auto& IY = property(4);
    const auto& J = property(5);
    const auto& IW = property(6);

    // uniform axial
    // strong axis bending near node
    // strong axis bending far node
    // weak axis bending near node
    // weak axis bending far node
    // torsion near node
    // torsion far node
    // warping near node
    // warping far node

    local_stiff.zeros(9, 9);
    local_stiff(0, 0) = E * A / length;

    auto factor = 2. * E * IZ / length;
    local_stiff(1, 2) = local_stiff(2, 1) = factor;
    factor *= 2.;
    local_stiff(1, 1) = local_stiff(2, 2) = factor;

    factor = 2. * E * IY / length;
    local_stiff(3, 4) = local_stiff(4, 3) = factor;
    factor *= 2.;
    local_stiff(3, 3) = local_stiff(4, 4) = factor;

    factor = (12. * E * IW / length / length + 1.2 * G * J) / length;
    local_stiff(5, 5) = local_stiff(6, 6) = factor;
    local_stiff(5, 6) = local_stiff(6, 5) = -factor;

    factor = 6. * E * IW / length / length + .1 * G * J;
    local_stiff(5, 7) = local_stiff(7, 5) = factor;
    local_stiff(5, 8) = local_stiff(8, 5) = factor;
    local_stiff(6, 7) = local_stiff(7, 6) = -factor;
    local_stiff(6, 8) = local_stiff(8, 6) = -factor;

    factor = 2. * E * IW / length - G * J * length / 30.;
    local_stiff(7, 8) = local_stiff(8, 7) = factor;

    factor = 4. * E * IW / length + 2. * G * J * length / 15.;
    local_stiff(7, 7) = local_stiff(8, 8) = factor;

    trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(local_stiff);

    ConstantStiffness(this);

    return SUANPAN_SUCCESS;
}

int EB31OS::update_status() {
    b_trans->update_status();

    const vec local_force = local_stiff * b_trans->to_local_vec(get_trial_displacement());

    trial_stiffness = b_trans->to_global_stiffness_mat(local_stiff);
    trial_resistance = b_trans->to_global_vec(local_force);

    if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(local_force);

    return SUANPAN_SUCCESS;
}

int EB31OS::commit_status() {
    b_trans->commit_status();
    return SUANPAN_SUCCESS;
}

int EB31OS::clear_status() {
    b_trans->clear_status();
    return SUANPAN_SUCCESS;
}

int EB31OS::reset_status() {
    b_trans->reset_status();
    return SUANPAN_SUCCESS;
}

std::vector<vec> EB31OS::record(const OutputType P) {
    if(P == OutputType::BEAME) return {b_trans->to_local_vec(get_current_displacement())};
    if(P == OutputType::BEAMS) return {vec{local_stiff * b_trans->to_local_vec(get_current_displacement())}};

    return {};
}

void EB31OS::print() {
    suanpan_info("A spatial beam element with warping DoF.\n");
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

vtkSmartPointer<vtkCell> EB31OS::Setup(const uvec& encoding) const {
    auto cell = vtkSmartPointer<vtkLine>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < b_node; ++I) {
        cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(encoding(I)));
        cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
    return cell;
}

mat EB31OS::GetData(const OutputType P) {
    if(OutputType::A == P) return resize(reshape(get_current_acceleration(), b_dof, b_node), 6, b_node);
    if(OutputType::V == P) return resize(reshape(get_current_velocity(), b_dof, b_node), 6, b_node);
    if(OutputType::U == P) return resize(reshape(get_current_displacement(), b_dof, b_node), 6, b_node);

    return {};
}

mat EB31OS::GetDeformation(const double amplifier) { return get_coordinate(3).t() + amplifier * reshape(get_current_displacement(), b_dof, b_node).eval().head_rows(3); }

#endif
