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

#include "Substepping.h"
#include <Domain/DomainBase.h>

Substepping::Substepping(const unsigned T, const unsigned MT, const unsigned MI)
    : Material(T)
    , max_iteration(MI)
    , mat_tag(MT) {}

Substepping::Substepping(const Substepping& old_obj)
    : Material(old_obj)
    , max_iteration(old_obj.max_iteration)
    , mat_tag(old_obj.mat_tag)
    , trial_mat_obj(suanpan::make_copy(old_obj.trial_mat_obj))
    , current_mat_obj(suanpan::make_copy(old_obj.current_mat_obj)) {}

int Substepping::initialize(const shared_ptr<DomainBase>& D) {
    current_mat_obj = suanpan::initialized_material_copy(D, mat_tag);

    if(nullptr == current_mat_obj) {
        suanpan_error("A valid host material is required.\n");
        return SUANPAN_FAIL;
    }

    trial_mat_obj = current_mat_obj->get_copy();

    PureWrapper(this);

    access::rw(density) = current_mat_obj->get_parameter(ParameterType::DENSITY);
    access::rw(material_type) = current_mat_obj->get_material_type();

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Substepping::get_copy() { return make_unique<Substepping>(*this); }

double Substepping::get_parameter(const ParameterType P) const { return current_mat_obj->get_parameter(P); }

const mat& Substepping::get_initial_damping() const { return current_mat_obj->get_initial_damping(); }

const vec& Substepping::get_initial_history() const { return current_mat_obj->get_initial_history(); }

const mat& Substepping::get_initial_stiffness() const { return current_mat_obj->get_initial_stiffness(); }

const mat& Substepping::get_current_damping() { return current_mat_obj->get_current_damping(); }

const mat& Substepping::get_current_secant() { return current_mat_obj->get_current_secant(); }

const mat& Substepping::get_current_stiffness() { return current_mat_obj->get_current_stiffness(); }

const mat& Substepping::get_trial_damping() { return trial_mat_obj->get_trial_damping(); }

const mat& Substepping::get_trial_secant() { return trial_mat_obj->get_trial_secant(); }

const mat& Substepping::get_trial_stiffness() { return trial_mat_obj->get_trial_stiffness(); }

const vec& Substepping::get_current_strain() { return current_mat_obj->get_current_strain(); }

const vec& Substepping::get_current_strain_rate() { return current_mat_obj->get_current_strain_rate(); }

const vec& Substepping::get_current_stress() { return current_mat_obj->get_current_stress(); }

const vec& Substepping::get_trial_strain() { return trial_mat_obj->get_trial_strain(); }

const vec& Substepping::get_trial_strain_rate() { return trial_mat_obj->get_trial_strain_rate(); }

const vec& Substepping::get_trial_stress() { return trial_mat_obj->get_trial_stress(); }

void Substepping::set_initial_history(const vec& H) {
    current_mat_obj->set_initial_history(H);
    trial_mat_obj->set_initial_history(H);
}

std::vector<vec> Substepping::record(const OutputType P) { return current_mat_obj->record(P); }

void Substepping::print() { if(current_mat_obj) current_mat_obj->print(); }

int Substepping::update_trial_status(const vec& t_strain) {
    incre_strain = t_strain - current_mat_obj->get_current_strain();

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_mat_obj = current_mat_obj->get_copy();

    auto accumulated_factor = 0., incre_factor = 1.;

    unsigned counter = 0, step = 0;
    while(true) {
        if(++counter == max_iteration) return SUANPAN_FAIL;

        if(SUANPAN_SUCCESS != trial_mat_obj->update_trial_status(current_mat_obj->get_current_strain() + (accumulated_factor + incre_factor) * incre_strain)) {
            step = 0;
            incre_factor *= .5;
            continue;
        }

        accumulated_factor += incre_factor;

        if(SUANPAN_SUCCESS != trial_mat_obj->commit_status()) return SUANPAN_FAIL;

        if(++step == 3) {
            step = 0;
            incre_factor *= 1.2;
        }

        if(fabs(incre_factor = std::min(incre_factor, 1. - accumulated_factor)) < tolerance) return SUANPAN_SUCCESS;
    }
}

int Substepping::clear_status() {
    current_mat_obj->clear_status();

    trial_mat_obj = current_mat_obj->get_copy();

    return SUANPAN_SUCCESS;
}

int Substepping::commit_status() {
    current_mat_obj = trial_mat_obj->get_copy();

    return SUANPAN_SUCCESS;
}

int Substepping::reset_status() {
    trial_mat_obj = current_mat_obj->get_copy();

    return SUANPAN_SUCCESS;
}
