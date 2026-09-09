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

#include "BoucWen.h"

#include <Toolbox/ridders.hpp>

BoucWen::BoucWen(const unsigned T, vec&& P)
    : DataBoucWen{.elastic_modulus = P(0), .yield_stress = P(1), .hardening = P(2), .beta = P(3), .n = P(4)}
    , Material1D(T, P(5)) {}

int BoucWen::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(1);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> BoucWen::unique_copy() { return std::make_unique<BoucWen>(*this); }

int BoucWen::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(std::fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    const auto n_strain = incre_strain(0) / yield_strain;

    trial_history = current_history;
    const auto& current_z = current_history(0); // z
    auto& z = trial_history(0);                 // z

    auto const_a{0.}, const_b{0.}, const_c{0.};
    if(std::fabs(n_strain) < .25) {
        // use Trapezoidal rule when step size is small
        const_a = (gamma + (current_z * n_strain >= 0. ? beta : -beta)) * std::pow(std::fabs(current_z), n);
        const_b = 2. * current_z + 2. * n_strain - n_strain * const_a;
        const_c = 2.;
    }
    else {
        // otherwise use backward Euler method
        const_a = 0.;
        const_b = current_z + n_strain;
        const_c = 1.;
    }

    auto counter = 0u;
    auto ref_error = 1.;
    auto try_bisection = false;
    while(true) {
        if(max_iteration == ++counter) {
            if(try_bisection) {
                suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
                return SUANPAN_FAIL;
            }

            try_bisection = true;
            counter = 2u;

            const auto approx_update = [&](const double in) {
                z = current_z + in;
                return const_c * z - const_b + n_strain * (gamma + (z * n_strain >= 0. ? beta : -beta)) * std::pow(std::max(datum::eps, std::fabs(z)), n);
            };

            ridders_guess(approx_update, 0., .25 * n_strain, tolerance);
        }

        const auto abs_z = std::max(datum::eps, std::fabs(z));
        const auto z_a = std::pow(abs_z, n - 1.);
        const auto b_term = gamma + (z * n_strain >= 0. ? beta : -beta);
        const auto p_term = b_term * z_a * abs_z;
        const auto t_term = n_strain * p_term;

        const auto residual = const_c * z - const_b + t_term;
        const auto jacobian = const_c + (z >= 0. ? n_strain : -n_strain) * b_term * n * z_a;
        const auto incre = residual / jacobian;
        const auto error = std::fabs(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);

        if(error < tolerance * ref_error || ((error < tolerance || std::fabs(residual) < tolerance) && counter > 5u)) {
            trial_stress = modulus_a * trial_strain + modulus_b * z;
            trial_stiffness = modulus_a + modulus_b / yield_strain * (const_c - p_term - const_a) / jacobian;

            return SUANPAN_SUCCESS;
        }

        z -= incre;
    }
}

void BoucWen::print() {
    suanpan_info("A Bouc-Wen material model.\n");
    Material1D::print();
}
