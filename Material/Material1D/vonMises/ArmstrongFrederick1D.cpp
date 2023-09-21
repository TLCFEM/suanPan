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

#include "ArmstrongFrederick1D.h"
#include <Recorder/OutputType.h>

ArmstrongFrederick1D::ArmstrongFrederick1D(const unsigned T, const double E, const double Y, const double S, const double H, const double M, vec&& A, vec&& B, const double R)
    : DataArmstrongFrederick1D{E, Y, S, H, M, std::forward<vec>(A), std::forward<vec>(B)}
    , Material1D(T, R) { access::rw(tolerance) = 1E-15; }

int ArmstrongFrederick1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(1 + size);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> ArmstrongFrederick1D::get_copy() { return make_unique<ArmstrongFrederick1D>(*this); }

double ArmstrongFrederick1D::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return elastic_modulus;
    return 0.;
}

int ArmstrongFrederick1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    auto& p = trial_history(size);

    auto yield_func = fabs(trial_stress(0) - accu(trial_history.head(size))) - std::max(0., yield + hardening * p + saturated * (1. - exp(-m * p)));

    if(yield_func < 0.) return SUANPAN_SUCCESS;

    auto gamma = 0.;
    double xi, jacobian;

    unsigned counter = 0;
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
        for(unsigned I = 0; I < size; ++I) {
            const auto denom = 1. + b(I) * gamma;
            sum_a += trial_history(I) / denom;
            sum_b += a(I) / denom;
        }

        yield_func = fabs(xi = trial_stress(0) - sum_a) - (elastic_modulus + sum_b) * gamma - k;

        jacobian = -elastic_modulus - dk;

        if(xi > 0.) for(unsigned I = 0; I < size; ++I) jacobian += (b(I) * trial_history(I) - a(I)) * pow(1. + b(I) * gamma, -2.);
        else for(unsigned I = 0; I < size; ++I) jacobian -= (b(I) * trial_history(I) + a(I)) * pow(1. + b(I) * gamma, -2.);

        const auto incre = yield_func / jacobian;
        suanpan_debug("Local iteration error: {:.5E}.\n", fabs(incre));
        if(fabs(incre) <= tolerance) break;

        gamma -= incre;
        p -= incre;
    }

    if(xi > 0.) {
        for(unsigned I = 0; I < size; ++I) trial_history(I) = (trial_history(I) + a(I) * gamma) / (1. + b(I) * gamma);

        trial_stress -= elastic_modulus * gamma;
    }
    else {
        for(unsigned I = 0; I < size; ++I) trial_history(I) = (trial_history(I) - a(I) * gamma) / (1. + b(I) * gamma);

        trial_stress += elastic_modulus * gamma;
    }

    trial_stiffness += elastic_modulus / jacobian * elastic_modulus;

    return SUANPAN_SUCCESS;
}

int ArmstrongFrederick1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int ArmstrongFrederick1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int ArmstrongFrederick1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

vector<vec> ArmstrongFrederick1D::record(const OutputType P) {
    if(P == OutputType::PEEQ) return {vec{current_history(size)}};

    return Material1D::record(P);
}

void ArmstrongFrederick1D::print() {
    suanpan_info("A uniaxial nonlinear hardening model using Armstrong-Frederick kinematic hardening rule.\n");
    Material1D::print();
}
