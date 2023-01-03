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

#include "SlipLock.h"
#include <Toolbox/utility.h>

SlipLock::SlipLock(const unsigned T, const double E, const double Y, const double H, const double R, const double D)
    : DataSlipLock{fabs(E), 1. / H, fabs(Y), fabs(R)}
    , Material1D(T, D) {}

int SlipLock::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> SlipLock::get_copy() { return make_unique<SlipLock>(*this); }

int SlipLock::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress;

    unsigned counter = 0;
    while(++counter < max_iteration) {
        const auto norm_stress = fabs(trial_stress(0) / yield_stress);
        const auto tmp_a = pow(norm_stress, R0) + 1.;
        const auto tmp_b = (1. - hardening_ratio) * pow(tmp_a, -1. / R0);
        trial_stiffness = elastic_modulus * tmp_a / (hardening_ratio * tmp_a + tmp_b);
        const auto error = trial_strain(0) - trial_stress(0) / elastic_modulus * (hardening_ratio + tmp_b);
        SP_D("Local iteration error: {:.5E}.\n", fabs(error));
        if(fabs(error) <= tolerance) return SUANPAN_SUCCESS;
        trial_stress += error * trial_stiffness;
        if(!suanpan::approx_equal(sign(trial_stress(0)), sign(trial_strain(0)))) trial_stress = 0.;
    }

    SP_E("Cannot converge within {} iterations.\n", max_iteration);
    return SUANPAN_FAIL;
}

int SlipLock::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_stiffness = initial_stiffness;
    return reset_status();
}

int SlipLock::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int SlipLock::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void SlipLock::print() {
    suanpan_info("A Slip-Lock material model using Menegotto-Pinto relationship.\n");
    Material1D::print();
}
