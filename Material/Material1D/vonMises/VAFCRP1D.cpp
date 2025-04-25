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

#include "VAFCRP1D.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

VAFCRP1D::VAFCRP1D(const unsigned T, DataVAFCRP1D&& D, const double R)
    : DataVAFCRP1D(std::move(D))
    , Material1D(T, R) { access::rw(tolerance) = 1E-15; }

int VAFCRP1D::initialize(const shared_ptr<DomainBase>& D) {
    if(nullptr != D) incre_time = &D->get_factory()->modify_incre_time();

    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(1 + size);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> VAFCRP1D::get_copy() { return make_unique<VAFCRP1D>(*this); }

double VAFCRP1D::get_parameter(const ParameterType P) const {
    if(ParameterType::ELASTICMODULUS == P) return elastic_modulus;
    return 0.;
}

int VAFCRP1D::update_trial_status(const vec& t_strain) {
    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * (incre_strain = (trial_strain = t_strain) - current_strain);

    trial_history = current_history;
    auto& p = trial_history(size);
    vec beta(&trial_history(0), size, false, true);

    if(fabs(trial_stress(0) - accu(beta)) < std::max(0., yield + hardening * p + saturated * (1. - exp(-m * p)))) return SUANPAN_SUCCESS;

    const auto norm_mu = mu / (incre_time && *incre_time > 0. ? *incre_time : 1.);

    auto gamma = 0., ref_error = 1.;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto exp_term = saturated * exp(-m * p);

        auto k = yield + saturated + hardening * p - exp_term;
        auto dk = hardening + m * exp_term;
        if(k < 0.) k = dk = 0.;

        const vec bottom = 1. + b * gamma;

        const auto xi = trial_stress(0) - accu(beta / bottom);
        const auto q = fabs(xi) - (elastic_modulus + accu(a / bottom)) * gamma;

        const auto fraction_term = norm_mu * gamma + 1.;
        const auto power_term = pow(fraction_term, epsilon - 1.);

        auto jacobian = -elastic_modulus - power_term * (norm_mu * epsilon * k + fraction_term * dk);

        if(xi > 0.) jacobian += accu((b % beta - a) / square(bottom));
        else jacobian -= accu((b % beta + a) / square(bottom));

        const auto residual = q - fraction_term * power_term * k;
        const auto incre = residual / jacobian;
        const auto error = fabs(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || fabs(residual) < tolerance * yield) && counter > 5u)) {
            if(xi > 0.) {
                beta += a * gamma;

                trial_stress -= elastic_modulus * gamma;
            }
            else {
                beta -= a * gamma;

                trial_stress += elastic_modulus * gamma;
            }
            beta /= bottom;

            trial_stiffness += elastic_modulus / jacobian * elastic_modulus;

            return SUANPAN_SUCCESS;
        }

        gamma -= incre;
        p -= incre;
    }
}

int VAFCRP1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int VAFCRP1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int VAFCRP1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void VAFCRP1D::print() {
    suanpan_info("A uniaxial VAFCRP material model.\n");
    Material1D::print();
}
