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

#include "LinearDamage.h"

#include <Toolbox/tensor.h>

const double LinearDamage::root_two_third = sqrt(2. / 3.);

LinearDamage::LinearDamage(const unsigned T, const unsigned MT, const double SE, const double EE, const double ED)
    : IsotropicDamage(T, MT)
    , s_strain(SE)
    , e_strain(EE)
    , e_damage(ED) {}

int LinearDamage::initialize(const shared_ptr<DomainBase>& D) {
    initialize_history(1);

    return IsotropicDamage::initialize(D);
}

void LinearDamage::compute_damage() {
    const auto dev_strain = tensor::dev(trial_strain);
    const auto eqv_strain = root_two_third * tensor::strain::norm(dev_strain);

    const auto compute_damage_factor = [&](const double t_strain) {
        if(t_strain <= s_strain) return vec{1., 0.};

        if(t_strain >= e_strain) return vec{e_damage, 0.};

        return vec{1. + slope * (t_strain - s_strain), slope};
    };

    trial_history = current_history;
    if(eqv_strain > trial_history(0)) {
        trial_history(0) = eqv_strain;
        const auto damage = compute_damage_factor(trial_history(0));
        const auto& d = damage(0);
        const auto& dd = damage(1);

        trial_stiffness *= d;

        trial_stiffness += trial_stress * (tensor::strain::norm_weight % dev_strain).t() * dd / eqv_strain / 1.5;

        trial_stress *= d;
    }
    else {
        const auto damage = compute_damage_factor(trial_history(0));
        const auto& d = damage(0);

        trial_stiffness *= d;

        trial_stress *= d;
    }
}

unique_ptr<Material> LinearDamage::get_copy() { return make_unique<LinearDamage>(*this); }
