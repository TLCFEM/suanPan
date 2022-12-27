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

#include "LinearDamage.h"
#include <Toolbox/tensorToolbox.h>

const double LinearDamage::root_two_third = sqrt(2. / 3.);

LinearDamage::LinearDamage(const unsigned T, const unsigned MT, const double SE, const double EE, const double ED)
    : IsotropicDamage(T, MT)
    , s_strain(SE)
    , e_strain(EE)
    , e_damage(ED) {}

void LinearDamage::compute_damage() {
    const auto dev_strain = tensor::dev(trial_strain);
    const auto eqv_strain = root_two_third * tensor::strain::norm(dev_strain);

    if(eqv_strain <= s_strain) return;

    if(eqv_strain >= e_strain) {
        trial_stiffness *= e_damage;
        trial_stress *= e_damage;
        return;
    }

    const auto damage = 1. + slope * (eqv_strain - s_strain);

    trial_stiffness *= damage;

    trial_stiffness += trial_stress * (tensor::strain::norm_weight % dev_strain).t() * slope / eqv_strain / 1.5;

    trial_stress *= damage;
}

unique_ptr<Material> LinearDamage::get_copy() { return make_unique<LinearDamage>(*this); }
