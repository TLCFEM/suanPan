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

#include "SteelBRB.h"

#include <Toolbox/ridders.hpp>
#include <Toolbox/utility.h>

pod2 SteelBRB::compute_t_yield_stress(const double plastic_strain) const {
    pod2 result;

    result[1] = t_const * std::exp(-plastic_strain / t_scalar);

    result[0] = t_saturated_stress - result[1] * t_scalar;

    return result;
}

pod2 SteelBRB::compute_c_yield_stress(const double plastic_strain) const {
    pod2 result;

    result[1] = c_const * std::exp(-plastic_strain / c_scalar);

    result[0] = c_saturated_stress - result[1] * c_scalar;

    return result;
}

SteelBRB::SteelBRB(const unsigned T, vec&& P)
    : DataSteelBRB{P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8)}
    , Material1D(T, P(9)) {}

int SteelBRB::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(2);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> SteelBRB::get_copy() { return std::make_unique<SteelBRB>(*this); }

int SteelBRB::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(std::fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = elastic_modulus) * incre_strain;

    trial_history = current_history;
    const auto& current_accumulated_strain = current_history(0); // u
    const auto& current_plastic_strain = current_history(1);     // \delta_1
    auto& accumulated_strain = trial_history(0);                 // u
    auto& plastic_strain = trial_history(1);                     // \delta_1

    // elastic unloading region
    if(std::signbit(trial_stress(0)) != std::signbit(incre_strain(0))) return SUANPAN_SUCCESS;

    const auto tension_flag = incre_strain(0) >= 0.;
    const auto& exponent = tension_flag ? t_exponent : c_exponent;
    const auto compute_stress = tension_flag ? std::mem_fn(&SteelBRB::compute_t_yield_stress) : std::mem_fn(&SteelBRB::compute_c_yield_stress);

    // need to identify two cases:
    // 1. same side full plastic loading
    // 2. unloading + reloading
    auto net_incre_strain = incre_strain(0);
    if(std::signbit(trial_stress(0)) != std::signbit(current_stress(0))) net_incre_strain += current_stress(0) / elastic_modulus;

    auto incre_plastic_strain = .5 * net_incre_strain;
    auto counter = 0u;
    auto ref_error = 1.;
    auto try_bisection = false;
    while(true) {
        if(max_iteration == ++counter) {
            const auto approx_update = [&](const auto in) {
                plastic_strain = current_plastic_strain + (incre_plastic_strain = in);
                trial_stress = elastic_modulus * (trial_strain - plastic_strain);

                const auto yield_stress = compute_stress(this, accumulated_strain = current_accumulated_strain + std::fabs(incre_plastic_strain))[0];
                return incre_plastic_strain - net_incre_strain * std::pow(std::fabs((trial_stress(0) - plastic_modulus * plastic_strain) / yield_stress), exponent);
            };

            const auto f1{approx_update(0.)}, f2{approx_update(net_incre_strain)};

            if(try_bisection || std::signbit(f1) == std::signbit(f2)) {
                suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
                return SUANPAN_FAIL;
            }

            try_bisection = true;
            counter = 2u;

            ridders(approx_update, 0., f1, net_incre_strain, f2, tolerance);
        }

        plastic_strain = current_plastic_strain + incre_plastic_strain;
        trial_stress = elastic_modulus * (trial_strain - plastic_strain);

        const auto sigma_y = compute_stress(this, accumulated_strain = current_accumulated_strain + std::fabs(incre_plastic_strain));
        const auto numerator = trial_stress(0) - plastic_modulus * plastic_strain;
        const auto fraction = numerator / sigma_y[0];
        const auto pow_term = std::pow(std::fabs(fraction), exponent);
        auto residual = -net_incre_strain * pow_term;
        const auto jacobian = 1. + exponent / numerator * residual * (s_modulus - fraction * (incre_plastic_strain >= 0. ? sigma_y[1] : -sigma_y[1]));
        residual += incre_plastic_strain;

        const auto incre = -residual / jacobian;
        const auto error = std::fabs(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);

        if(error < tolerance * ref_error || ((error < tolerance || std::fabs(residual) < tolerance) && counter > 5u)) {
            trial_stiffness *= 1. - (pow_term + net_incre_strain * elastic_modulus * exponent * pow_term / numerator) / jacobian;

            return SUANPAN_SUCCESS;
        }

        incre_plastic_strain = suanpan::clamp(incre_plastic_strain + incre, net_incre_strain, 0.);
    }
}

int SteelBRB::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int SteelBRB::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    current_history = trial_history;
    return SUANPAN_SUCCESS;
}

int SteelBRB::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    trial_history = current_history;
    return SUANPAN_SUCCESS;
}

void SteelBRB::print() {
    suanpan_info("A steel model for BRB. doi:10.1016/j.jcsr.2011.07.017\n");
    Material1D::print();
}
