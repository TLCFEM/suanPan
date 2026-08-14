/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#include "CoupleMaterial.h"

#include <Domain/DomainBase.h>

void ConstantStiffness(CoupleMaterial* M) {
    M->current_stiffness = mat(M->initial_stiffness.memptr(), M->initial_stiffness.n_rows, M->initial_stiffness.n_cols, false, true);
    M->trial_stiffness = mat(M->initial_stiffness.memptr(), M->initial_stiffness.n_rows, M->initial_stiffness.n_cols, false, true);
}

CoupleMaterial::CoupleMaterial(const double L)
    : DataCoupleMaterial{.characteristic_length = L} {}

const vec& CoupleMaterial::get_trial_curvature() { return trial_curvature; }

const vec& CoupleMaterial::get_trial_moment() { return trial_moment; }

const mat& CoupleMaterial::get_trial_stiffness() { return trial_stiffness; }

const vec& CoupleMaterial::get_current_curvature() { return current_curvature; }

const vec& CoupleMaterial::get_current_moment() { return current_moment; }

const mat& CoupleMaterial::get_current_stiffness() { return current_stiffness; }

const mat& CoupleMaterial::get_initial_stiffness() const { return initial_stiffness; }

int CoupleMaterial::update_incre_status(const double i_strain) { return update_incre_status(vec{i_strain}); }

int CoupleMaterial::update_incre_status(const double i_strain, const double i_strain_rate) { return update_incre_status(vec{i_strain}, vec{i_strain_rate}); }

int CoupleMaterial::update_incre_status(const double i_strain, const double i_strain_rate, const double i_strain_acc) { return update_incre_status(vec{i_strain}, vec{i_strain_rate}, vec{i_strain_acc}); }

int CoupleMaterial::update_trial_status(const double t_strain) { return update_trial_status(vec{t_strain}); }

int CoupleMaterial::update_trial_status(const double t_strain, const double t_strain_rate) { return update_trial_status(vec{t_strain}, vec{t_strain_rate}); }

int CoupleMaterial::update_trial_status(const double t_strain, const double t_strain_rate, const double t_strain_acc) { return update_trial_status(vec{t_strain}, vec{t_strain_rate}, vec{t_strain_acc}); }

int CoupleMaterial::update_incre_status(const vec& i_curvature) { return update_trial_status(current_curvature + i_curvature); }

int CoupleMaterial::update_incre_status(const vec& i_curvature, const vec&) { return update_trial_status(current_curvature + i_curvature); }

int CoupleMaterial::update_incre_status(const vec& i_curvature, const vec&, const vec&) { return update_trial_status(current_curvature + i_curvature); }

int CoupleMaterial::update_trial_status(const vec& t_curvature) {
    trial_moment = trial_stiffness * (trial_curvature = t_curvature);
    return SUANPAN_SUCCESS;
}

int CoupleMaterial::update_trial_status(const vec& t_curvature, const vec&) { return update_trial_status(t_curvature); }

int CoupleMaterial::update_trial_status(const vec& t_curvature, const vec&, const vec&) { return update_trial_status(t_curvature); }

int CoupleMaterial::clear_status() {
    if(!current_curvature.is_empty()) current_curvature.zeros();
    if(!current_moment.is_empty()) current_moment.zeros();

    if(!trial_curvature.is_empty()) trial_curvature.zeros();
    if(!trial_moment.is_empty()) trial_moment.zeros();

    if(!initial_stiffness.is_empty()) trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int CoupleMaterial::commit_status() {
    if(!trial_curvature.is_empty()) current_curvature = trial_curvature;
    if(!trial_moment.is_empty()) current_moment = trial_moment;
    if(!trial_stiffness.is_empty()) current_stiffness = trial_stiffness;

    return SUANPAN_SUCCESS;
}

int CoupleMaterial::reset_status() {
    if(!trial_curvature.is_empty()) trial_curvature = current_curvature;
    if(!trial_moment.is_empty()) trial_moment = current_moment;
    if(!trial_stiffness.is_empty()) trial_stiffness = current_stiffness;

    return SUANPAN_SUCCESS;
}

IsotropicCouple::IsotropicCouple(const double E, const double V, const double L)
    : CoupleMaterial(L)
    , elastic_modulus(E)
    , poissons_ratio(V) {}

void IsotropicCouple::initialize(const shared_ptr<DomainBase>&) {
    if(characteristic_length < 0.) {
        access::rw(characteristic_length) = 1.;
        suanpan_warning("Characteristic length is set to unity.\n");
    }

    initial_stiffness = 2. * characteristic_length * characteristic_length * elastic_modulus / (1. + poissons_ratio) * eye(2, 2);

    trial_curvature = current_curvature.zeros(2);
    trial_moment = current_moment.zeros(2);

    ConstantStiffness(this);
}
