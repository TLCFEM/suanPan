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

#include "BoucWen.h"

BoucWen::BoucWen(const unsigned T, vec&& P)
    : DataBoucWen{P(0), P(1), P(2), P(3), P(4)}
    , Material1D(T, P(5)) {}

int BoucWen::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(1);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> BoucWen::get_copy() { return make_unique<BoucWen>(*this); }

int BoucWen::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    const auto n_strain = incre_strain(0) / yield_strain;

    trial_history = current_history;
    const auto& current_z = current_history(0); // z
    auto& z = trial_history(0);                 // z

    auto incre = .5 * n_strain;
    unsigned counter = 0;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("BoucWen cannot converge within %u iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        z += incre;

        const auto p_term = (gamma + (z * n_strain >= 0. ? beta : -beta)) * pow(std::max(datum::eps, fabs(z)), n);
        const auto t_term = n_strain * p_term;

        const auto residual = z - current_z + t_term - n_strain;
        const auto jacobian = z + n * t_term;

        const auto error = fabs(incre = -residual * z / jacobian);

        suanpan_debug("BoucWen local iteration error: %.5E.\n", error);

        if(error <= tolerance) {
            trial_stress = modulus_a * trial_strain + modulus_b * z;
            trial_stiffness = modulus_a + modulus_b / yield_strain * (1. - p_term) * z / jacobian;

            return SUANPAN_SUCCESS;
        }
    }
}

int BoucWen::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int BoucWen::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    current_history = trial_history;
    return SUANPAN_SUCCESS;
}

int BoucWen::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    trial_history = current_history;
    return SUANPAN_SUCCESS;
}

void BoucWen::print() {
    suanpan_info("Bouc-Wen model.\n");
    Material1D::print();
}
