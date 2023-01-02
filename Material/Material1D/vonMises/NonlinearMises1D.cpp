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

#include "NonlinearMises1D.h"

NonlinearMises1D::NonlinearMises1D(const unsigned T, const double E, const double R)
    : DataMises1D{fabs(E)}
    , Material1D(T, R) {}

int NonlinearMises1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(2);

    return SUANPAN_SUCCESS;
}

int NonlinearMises1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& plastic_strain = trial_history(0);
    auto& back_stress = trial_history(1);

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    const auto shifted_stress = trial_stress(0) - back_stress;
    const auto norm_shifted_stress = fabs(shifted_stress);

    if(auto yield_func = norm_shifted_stress - std::max(0., compute_k(plastic_strain)); yield_func >= 0.) {
        const auto current_h = compute_h(plastic_strain);
        auto gamma = 0., incre_h = 0.;
        double dkdh;
        auto counter = 0u;
        while(true) {
            if(max_iteration == ++counter) {
                SP_E("Cannot converge within {} iterations.\n", max_iteration);
                return SUANPAN_FAIL;
            }

            const auto incre_gamma = yield_func / (elastic_modulus + (dkdh = compute_dk(plastic_strain) + compute_dh(plastic_strain)));
            const auto abs_error = fabs(incre_gamma);
            suanpan_extra_debug("NonlinearMises1D local iteration error: %.5E.\n", abs_error);
            if(abs_error <= tolerance) break;
            incre_h = compute_h(plastic_strain = current_history(0) + (gamma += incre_gamma)) - current_h;
            yield_func = norm_shifted_stress - elastic_modulus * gamma - std::max(0., compute_k(plastic_strain)) - incre_h;
        }

        if(shifted_stress > 0.) {
            back_stress += incre_h;
            trial_stress -= elastic_modulus * gamma;
        }
        else {
            back_stress -= incre_h;
            trial_stress += elastic_modulus * gamma;
        }

        trial_stiffness *= dkdh / (dkdh + elastic_modulus);
    }

    return SUANPAN_SUCCESS;
}

int NonlinearMises1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int NonlinearMises1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int NonlinearMises1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}
