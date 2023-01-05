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

#include "Rotation2D.h"
#include <Domain/DomainBase.h>
#include <Toolbox/tensor.h>

Rotation2D::Rotation2D(const unsigned T, const unsigned MT, const double A)
    : Material2D(T, PlaneType::N, 0.)
    , mat_tag(MT)
    , trans_mat(transform::strain::trans(A)) {}

Rotation2D::Rotation2D(const Rotation2D& old_obj)
    : Material2D(old_obj)
    , mat_tag(old_obj.mat_tag)
    , mat_obj(suanpan::make_copy(old_obj.mat_obj))
    , trans_mat(old_obj.trans_mat) {}

int Rotation2D::initialize(const shared_ptr<DomainBase>& D) {
    mat_obj = suanpan::initialized_material_copy(D, mat_tag);

    if(nullptr == mat_obj || mat_obj->get_material_type() != MaterialType::D2) {
        suanpan_error("A valid 2D host material is required.\n");
        return SUANPAN_FAIL;
    }

    access::rw(density) = mat_obj->get_parameter(ParameterType::DENSITY);

    trial_stiffness = current_stiffness = initial_stiffness = trans_mat.t() * mat_obj->get_initial_stiffness() * trans_mat;

    return SUANPAN_SUCCESS;
}

double Rotation2D::get_parameter(const ParameterType P) const { return mat_obj->get_parameter(P); }

unique_ptr<Material> Rotation2D::get_copy() { return make_unique<Rotation2D>(*this); }

int Rotation2D::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;

    if(SUANPAN_SUCCESS != mat_obj->update_trial_status(trans_mat * trial_strain)) return SUANPAN_FAIL;

    trial_stress = trans_mat.t() * mat_obj->get_trial_stress();
    trial_stiffness = trans_mat.t() * mat_obj->get_trial_stiffness() * trans_mat;

    return SUANPAN_SUCCESS;
}

int Rotation2D::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    return mat_obj->clear_status();
}

int Rotation2D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return mat_obj->commit_status();
}

int Rotation2D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return mat_obj->reset_status();
}

vector<vec> Rotation2D::record(const OutputType P) { return mat_obj->record(P); }

void Rotation2D::print() {
    suanpan_info("A rotation wrapper with the underlying material.\n");
    if(nullptr != mat_obj) mat_obj->print();
}
