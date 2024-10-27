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

#include "ArmstrongFrederick1D.h"

ArmstrongFrederick1D::ArmstrongFrederick1D(const unsigned T, DataArmstrongFrederick1D&& D, const double R)
    : DataArmstrongFrederick1D(std::move(D))
    , Material1D(T, R) {}

int ArmstrongFrederick1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(size + 4);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> ArmstrongFrederick1D::get_copy() { return make_unique<ArmstrongFrederick1D>(*this); }

double ArmstrongFrederick1D::get_parameter(const ParameterType P) const {
    if(ParameterType::ELASTICMODULUS == P) return elastic_modulus;
    return 0.;
}

int ArmstrongFrederick1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    auto& q = trial_history(size);
    const auto& current_r = current_history(size + 2);
    const auto& current_theta = current_history(size + 3);
    auto& ep = trial_history(size + 1);
    auto& r = trial_history(size + 2);
    auto& theta = trial_history(size + 3);

    auto yield_func = fabs(trial_stress(0) - accu(trial_history.head(size))) - std::max(0., yield + hardening * q + saturation * (1. - exp(-ms * q)) - reduction * (1. - exp(-mr * r)));

    if(yield_func < 0.) return SUANPAN_SUCCESS;

    auto gamma = 0.;
    double xi, jacobian, dr = 0.;

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto s_term = saturation * exp(-ms * q);
        const auto r_term = reduction * exp(-mr * r);

        auto k = yield + saturation - reduction + hardening * q - s_term + r_term;
        auto dk = hardening + ms * s_term - mr * r_term * dr;
        if(k < 0.) k = dk = 0.;

        auto sum_a = 0., sum_b = 0.;
        for(auto I = 0u; I < size; ++I) {
            const auto denom = 1. + b(I) * gamma;
            sum_a += trial_history(I) / denom;
            sum_b += a(I) / denom;
        }

        yield_func = fabs(xi = trial_stress(0) - sum_a) - (elastic_modulus + sum_b) * gamma - k;

        jacobian = -elastic_modulus - dk;

        if(xi > 0.) for(auto I = 0u; I < size; ++I) jacobian += (b(I) * trial_history(I) - a(I)) * pow(1. + b(I) * gamma, -2.);
        else for(auto I = 0u; I < size; ++I) jacobian -= (b(I) * trial_history(I) + a(I)) * pow(1. + b(I) * gamma, -2.);

        const auto incre = yield_func / jacobian;
        const auto error = fabs(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || fabs(yield_func) < tolerance) && counter > 5u)) break;

        gamma -= incre;
        q -= incre;
        ep -= xi > 0. ? incre : -incre;

        r = current_r;
        theta = current_theta;

        if(const auto h = fabs(ep - current_theta) - current_r; h > 0.) {
            const auto nh = ep > current_theta ? 1. : -1.;
            r += memory * h;
            theta += nh * (1. - memory) * h;
            dr = xi > 0. ? memory * nh : -memory * nh;
        }
        else dr = 0.;
    }

    if(xi > 0.) {
        for(auto I = 0u; I < size; ++I) trial_history(I) = (trial_history(I) + a(I) * gamma) / (1. + b(I) * gamma);

        trial_stress -= elastic_modulus * gamma;
    }
    else {
        for(auto I = 0u; I < size; ++I) trial_history(I) = (trial_history(I) - a(I) * gamma) / (1. + b(I) * gamma);

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

void ArmstrongFrederick1D::print() {
    suanpan_info("A uniaxial nonlinear hardening model using Armstrong-Frederick kinematic hardening rule.\n");
    Material1D::print();
}
