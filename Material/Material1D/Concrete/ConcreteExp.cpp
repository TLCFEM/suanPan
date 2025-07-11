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

#include "ConcreteExp.h"

pod2 ConcreteExp::compute_compression_initial_reverse() const {
    pod2 response;

    response[1] = middle_point * f_c;
    response[0] = response[1] / elastic_modulus;

    return response;
}

pod2 ConcreteExp::compute_tension_initial_reverse() const {
    pod2 response;

    response[1] = middle_point * f_t;
    response[0] = response[1] / elastic_modulus;

    return response;
}

pod2 ConcreteExp::compute_compression_backbone(const double n_strain) const {
    pod2 response;

    response[0] = (response[1] = elastic_modulus) * n_strain;

    if(n_strain * elastic_modulus >= f_c) return response;

    double exp_term, jacobian;
    auto stress = .25 * f_c / a_c * pow(1. + a_c, 2.);

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return response;
        }

        exp_term = std::exp(b_c * (n_strain - stress / elastic_modulus));
        const auto residual = (1. + a_c - a_c * exp_term) * exp_term - stress / f_c;
        jacobian = b_c / elastic_modulus * (2. * a_c * exp_term - 1. - a_c) * exp_term - 1. / f_c;
        const auto incre = residual / jacobian;
        const auto error = std::fabs(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local compression iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || std::fabs(residual) < tolerance) && counter > 5u)) break;

        stress -= incre;
    }

    response[0] = stress;
    response[1] = b_c * (2. * a_c * exp_term - 1. - a_c) * exp_term / jacobian;

    return response;
}

pod2 ConcreteExp::compute_tension_backbone(const double n_strain) const {
    pod2 response;

    response[0] = (response[1] = elastic_modulus) * n_strain;

    if(n_strain * elastic_modulus <= f_t) return response;

    double exp_term, jacobian;
    auto stress = f_t;

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return response;
        }

        exp_term = std::exp(-b_t * (n_strain - stress / elastic_modulus));
        const auto residual = (1. + a_t - a_t * exp_term) * exp_term - stress / f_t;
        jacobian = b_t / elastic_modulus * (1. + a_t - 2. * a_t * exp_term) * exp_term - 1. / f_t;
        const auto incre = residual / jacobian;
        const auto error = std::fabs(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local tension iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || std::fabs(residual) < tolerance) && counter > 5u)) break;

        stress -= incre;
    }

    response[0] = stress;
    response[1] = b_t * (1. + a_t - 2. * a_t * exp_term) * exp_term / jacobian;

    return response;
}

double ConcreteExp::compute_compression_residual(const double reverse_c_strain, const double reverse_c_stress) const { return std::min(0., reverse_c_strain - reverse_c_stress * (reverse_c_strain / f_c + .57 / elastic_modulus) / (reverse_c_stress / f_c + .57)); }

double ConcreteExp::compute_tension_residual(const double reverse_t_strain, const double reverse_t_stress) const { return std::max(0., reverse_t_strain - reverse_t_stress * (reverse_t_strain / f_t + .67 / elastic_modulus) / (reverse_t_stress / f_t + .67)); }

ConcreteExp::ConcreteExp(const unsigned T, const double E, const double FT, const double AT, const double GT, const double FC, const double AC, const double GC, const double M, const double R)
    : DataConcreteExp{E, std::fabs(FT), -std::fabs(FC) * 4. * AC * std::pow(1. + AC, -2.), AT, AC, std::fabs(FT) / GT * (1. + .5 * AT), std::fabs(FC) * 4. * AC * std::pow(1. + AC, -2.) / GC * (1. + .5 * AC)}
    , SimpleHysteresis(T, M, R) {}

int ConcreteExp::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    return SUANPAN_SUCCESS;
}

double ConcreteExp::get_parameter(const ParameterType P) const {
    if(ParameterType::ELASTICMODULUS == P) return elastic_modulus;
    if(ParameterType::PEAKSTRAIN == P) return f_c / elastic_modulus;
    if(ParameterType::CRACKSTRAIN == P) return f_t / elastic_modulus;
    return 0.;
}

unique_ptr<Material> ConcreteExp::get_copy() { return std::make_unique<ConcreteExp>(*this); }
