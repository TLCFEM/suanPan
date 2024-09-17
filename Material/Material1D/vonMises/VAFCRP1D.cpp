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

#include "VAFCRP1D.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Recorder/OutputType.h>

constexpr double VAFCRP1D::unit_time = 1.;

VAFCRP1D::VAFCRP1D(const unsigned T, DataVAFCRP1D&& D, const double R)
    : DataVAFCRP1D(std::move(D))
    , Material1D(T, R) { access::rw(tolerance) = 1E-15; }

int VAFCRP1D::initialize(const shared_ptr<DomainBase>& D) {
    incre_time = D == nullptr ? &unit_time : &D->get_factory()->modify_incre_time();

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

    if(fabs(trial_stress(0) - accu(trial_history.head(size))) < std::max(0., yield + hardening * p + saturated * (1. - exp(-m * p)))) return SUANPAN_SUCCESS;

    auto gamma = 0.;
    double xi, jacobian, exp_gamma;

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto exp_term = saturated * exp(-m * p);

        auto k = yield + saturated + hardening * p - exp_term;
        auto dk = hardening + m * exp_term;
        if(k < 0.) k = dk = 0.;

        auto sum_a = 0., sum_b = 0.;
        for(auto I = 0u; I < size; ++I) {
            const auto denom = 1. + b(I) * gamma;
            sum_a += trial_history(I) / denom;
            sum_b += a(I) / denom;
        }

        const auto q = fabs(xi = trial_stress(0) - sum_a) - (elastic_modulus + sum_b) * gamma;

        exp_gamma = pow(*incre_time / (*incre_time + mu * gamma), epsilon);

        jacobian = -elastic_modulus - epsilon * mu * q / (*incre_time + mu * gamma);

        if(xi > 0.) for(auto I = 0u; I < size; ++I) jacobian += (b(I) * trial_history(I) - a(I)) * pow(1. + b(I) * gamma, -2.);
        else for(auto I = 0u; I < size; ++I) jacobian -= (b(I) * trial_history(I) + a(I)) * pow(1. + b(I) * gamma, -2.);

        const auto residual = q * exp_gamma - k;
        const auto incre = residual / ((jacobian *= exp_gamma) -= dk);
        const auto error = fabs(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || (fabs(residual) < tolerance && counter > 5u)) break;

        gamma -= incre;
        p -= incre;
    }

    if(xi > 0.) {
        for(auto I = 0u; I < size; ++I) trial_history(I) = (trial_history(I) + a(I) * gamma) / (1. + b(I) * gamma);

        trial_stress -= elastic_modulus * gamma;
    }
    else {
        for(auto I = 0u; I < size; ++I) trial_history(I) = (trial_history(I) - a(I) * gamma) / (1. + b(I) * gamma);

        trial_stress += elastic_modulus * gamma;
    }

    trial_stiffness += elastic_modulus / jacobian * elastic_modulus * exp_gamma;

    return SUANPAN_SUCCESS;
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
