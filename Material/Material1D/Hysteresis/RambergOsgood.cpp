/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "RambergOsgood.h"
#include <Toolbox/utility.h>

RambergOsgood::RambergOsgood(const unsigned T, const double E, const double Y, const double O, const double N, const double R)
    : DataRambergOsgood{fabs(E), fabs(Y), fabs(O), fabs(N)}
    , Material1D(T, R) {}

int RambergOsgood::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(5);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> RambergOsgood::get_copy() { return make_unique<RambergOsgood>(*this); }

int RambergOsgood::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& load_sign = trial_history(0);
    auto& reverse_strain = trial_history(1);
    auto& reverse_stress = trial_history(2);
    auto& p_reverse_strain = trial_history(3);
    auto& p_reverse_stress = trial_history(4);

    if(const auto trial_load_sign = suanpan::sign(incre_strain(0)); !suanpan::approx_equal(trial_load_sign, load_sign)) {
        if(!suanpan::approx_equal(load_sign, 0.)) {
            p_reverse_strain = reverse_strain;
            p_reverse_stress = reverse_stress;
            reverse_strain = current_strain(0);
            reverse_stress = current_stress(0);
        }
        load_sign = trial_load_sign;
    }

    const auto elastic_predictor = elastic_modulus * fabs(trial_strain(0) - reverse_strain);

    const auto norm_yield_stress = std::max(datum::eps, 0. == p_reverse_stress || fabs(p_reverse_stress) < fabs(reverse_stress) ? std::max(fabs(reverse_stress - p_reverse_stress), yield_stress + fabs(reverse_stress)) : fabs(reverse_stress - p_reverse_stress));

    const auto pow_a = pow(norm_yield_stress, nm);

    auto norm_stress = fabs(current_stress(0) - reverse_stress);

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto pow_b = offset * pow(norm_stress, nm);
        const auto residual = norm_stress * (pow_a + pow_b) - elastic_predictor * pow_a;
        const auto jacobian = pow_a + n * pow_b;
        const auto incre = residual / jacobian;

        const auto error = fabs(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);

        if(error < tolerance * ref_error || (fabs(residual) < tolerance && counter > 5u)) {
            trial_stress = load_sign * norm_stress + reverse_stress;
            trial_stiffness = elastic_modulus * pow_a / jacobian;
            return SUANPAN_SUCCESS;
        }
        norm_stress -= incre;
    }
}

int RambergOsgood::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int RambergOsgood::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int RambergOsgood::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void RambergOsgood::print() {
    suanpan_info("A Ramberg-Osgood material model.\n");
    Material1D::print();
}
