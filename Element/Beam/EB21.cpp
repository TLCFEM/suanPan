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

#include "EB21.h"

#include <Domain/DomainBase.h>
#include <Element/Utility/B2DC.h>
#include <Material/Material1D/Material1D.h>
#include <Recorder/OutputType.h>

EB21::EB21(const unsigned T, uvec&& N, const double A, const double I, const unsigned M, const bool F)
    : MaterialElement1D(T, b_node, b_dof, std::move(N), uvec{M}, F, {Node::DOF::U1, Node::DOF::U2, Node::DOF::UR3})
    , area(A)
    , moment_inertia(I)
    , b_trans(F ? std::make_unique<B2DC>() : std::make_unique<B2DL>()) {}

int EB21::initialize(const shared_ptr<DomainBase>& D) {
    b_trans->set_element_ptr(this);

    b_material = D->get<Material>(material_tag(0))->get_copy();

    // stiffness
    const auto tmp_a = as_scalar(b_material->get_initial_stiffness()) / b_trans->get_length();

    local_stiff.zeros(3, 3);
    local_stiff(0, 0) = tmp_a * area;
    local_stiff(1, 1) = local_stiff(2, 2) = 2. * (local_stiff(1, 2) = local_stiff(2, 1) = 2. * tmp_a * moment_inertia);

    trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(local_stiff);

    if(b_material->get_density() > 0.) initial_mass = b_trans->to_global_mass_mat(b_material->get_density() * area);

    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

int EB21::update_status() {
    b_trans->update_status();

    const vec local_force = local_stiff * b_trans->to_local_vec(get_trial_displacement());

    trial_stiffness = b_trans->to_global_stiffness_mat(local_stiff);
    trial_resistance = b_trans->to_global_vec(local_force);

    if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(local_force);

    return SUANPAN_SUCCESS;
}

int EB21::commit_status() {
    b_trans->commit_status();
    return b_material->commit_status();
}

int EB21::clear_status() {
    b_trans->clear_status();
    return b_material->clear_status();
}

int EB21::reset_status() {
    b_trans->reset_status();
    return b_material->reset_status();
}

std::vector<vec> EB21::record(const OutputType P) const {
    if(P == OutputType::BEAME) return {b_trans->to_local_vec(get_current_displacement())};
    if(P == OutputType::BEAMS) return {vec{local_stiff * b_trans->to_local_vec(get_current_displacement())}};

    return {};
}

void EB21::print() {
    suanpan_info("An elastic B21 element{}", nlgeom ? " with corotational formulation.\n" : ".\n");
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

vtkSmartPointer<vtkCell> EB21::GetCell() const { return vtkSmartPointer<vtkLine>::New(); }

mat EB21::GetData(const OutputType P) {
    const auto remap = [&](vec&& in) {
        mat data(6, b_node, fill::zeros);
        data.rows(uvec{0, 1, 5}) = reshape(in, b_dof, b_node);
        return data;
    };

    if(OutputType::A == P) return remap(get_current_acceleration());
    if(OutputType::V == P) return remap(get_current_velocity());
    if(OutputType::U == P) return remap(get_current_displacement());

    return {};
}

mat EB21::GetDeformation(const double amplifier) { return get_coordinate(2).t() + amplifier * reshape(get_current_displacement(), b_dof, b_node).eval().head_rows(2); }

#endif
