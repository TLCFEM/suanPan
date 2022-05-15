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

#include "CoulombFriction.h"

CoulombFriction::CoulombFriction(const unsigned T, const double F, const double S)
    : DataCoulombFriction{2. / datum::pi * fabs(F), fabs(S)}
    , Material1D(T, 0.) {}

int CoulombFriction::initialize(const shared_ptr<DomainBase>&) {
    trial_damping = current_damping = initial_damping = friction_force * factor;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> CoulombFriction::get_copy() { return make_unique<CoulombFriction>(*this); }

int CoulombFriction::update_trial_status(const vec&) {
    suanpan_error("CoulombFriction receives strain only from the associated element, check the model.\n");
    return SUANPAN_FAIL;
}

int CoulombFriction::update_trial_status(const vec&, const vec& t_strain_rate) {
    incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

    if(fabs(incre_strain_rate(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = friction_force * atan(factor * trial_strain_rate(0));

    trial_damping = friction_force * factor / (1. + pow(factor * trial_strain_rate(0), 2.));

    return SUANPAN_SUCCESS;
}

int CoulombFriction::clear_status() {
    current_strain_rate.zeros();
    current_stress.zeros();
    current_damping = initial_damping;
    return reset_status();
}

int CoulombFriction::commit_status() {
    current_strain_rate = trial_strain_rate;
    current_stress = trial_stress;
    current_damping = trial_damping;
    return SUANPAN_SUCCESS;
}

int CoulombFriction::reset_status() {
    trial_strain_rate = current_strain_rate;
    trial_stress = current_stress;
    trial_damping = current_damping;
    return SUANPAN_SUCCESS;
}

void CoulombFriction::print() { suanpan_info("A Coulomb friction model.\n"); }
