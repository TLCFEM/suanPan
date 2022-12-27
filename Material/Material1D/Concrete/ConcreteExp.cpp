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

#include "ConcreteExp.h"

podarray<double> ConcreteExp::compute_compression_initial_reverse() const {
    podarray<double> response(2);

    response(1) = middle_point * f_c;
    response(0) = response(1) / elastic_modulus;

    return response;
}

podarray<double> ConcreteExp::compute_tension_initial_reverse() const {
    podarray<double> response(2);

    response(1) = middle_point * f_t;
    response(0) = response(1) / elastic_modulus;

    return response;
}

podarray<double> ConcreteExp::compute_compression_backbone(const double n_strain) const {
    podarray<double> response(2);

    response(0) = (response(1) = elastic_modulus) * n_strain;

    if(n_strain * elastic_modulus >= f_c) return response;

    auto counter = 0;

    double exp_term, jacobian;
    auto stress = .25 * f_c / a_c * pow(1. + a_c, 2.);

    while(true) {
        if(++counter == max_iteration) {
            suanpan_error("ConcreteExp cannot converge within %u iterations.\n", max_iteration);
            return response;
        }

        exp_term = exp(b_c * (n_strain - stress / elastic_modulus));
        const auto residual = (1. + a_c - a_c * exp_term) * exp_term - stress / f_c;
        jacobian = b_c / elastic_modulus * (2. * a_c * exp_term - 1. - a_c) * exp_term - 1. / f_c;
        const auto incre = residual / jacobian;

        suanpan_debug("ConcreteExp local compression iteration error: %.5E.\n", fabs(incre));
        if(fabs(incre) < tolerance) break;

        stress -= incre;
    }

    response(0) = stress;
    response(1) = b_c * (2. * a_c * exp_term - 1. - a_c) * exp_term / jacobian;

    return response;
}

podarray<double> ConcreteExp::compute_tension_backbone(const double n_strain) const {
    podarray<double> response(2);

    response(0) = (response(1) = elastic_modulus) * n_strain;

    if(n_strain * elastic_modulus <= f_t) return response;

    auto counter = 0;

    double exp_term, jacobian;
    auto stress = f_t;

    while(true) {
        if(++counter == max_iteration) {
            suanpan_error("ConcreteExp cannot converge within %u iterations.\n", max_iteration);
            return response;
        }

        exp_term = exp(-b_t * (n_strain - stress / elastic_modulus));
        const auto residual = (1. + a_t - a_t * exp_term) * exp_term - stress / f_t;
        jacobian = b_t / elastic_modulus * (1. + a_t - 2. * a_t * exp_term) * exp_term - 1. / f_t;
        const auto incre = residual / jacobian;

        suanpan_debug("ConcreteExp local tension iteration error: %.5E.\n", fabs(incre));
        if(fabs(incre) < tolerance) break;

        stress -= incre;
    }

    response(0) = stress;
    response(1) = b_t * (1. + a_t - 2. * a_t * exp_term) * exp_term / jacobian;

    return response;
}

double ConcreteExp::compute_compression_residual(const double reverse_c_strain, const double reverse_c_stress) const { return std::min(0., reverse_c_strain - reverse_c_stress * (reverse_c_strain / f_c + .57 / elastic_modulus) / (reverse_c_stress / f_c + .57)); }

double ConcreteExp::compute_tension_residual(const double reverse_t_strain, const double reverse_t_stress) const { return std::max(0., reverse_t_strain - reverse_t_stress * (reverse_t_strain / f_t + .67 / elastic_modulus) / (reverse_t_stress / f_t + .67)); }

ConcreteExp::ConcreteExp(const unsigned T, const double E, const double FT, const double AT, const double GT, const double FC, const double AC, const double GC, const double M, const double R)
    : DataConcreteExp{E, fabs(FT), -fabs(FC) * 4. * AC * pow(1. + AC, -2.), AT, AC, fabs(FT) / GT * (1. + .5 * AT), fabs(FC) * 4. * AC * pow(1. + AC, -2.) / GC * (1. + .5 * AC)}
    , SimpleHysteresis(T, M, R) { access::rw(tolerance) = 1E-13; }

int ConcreteExp::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    return SUANPAN_SUCCESS;
}

double ConcreteExp::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return elastic_modulus;
    if(ParameterType::PEAKSTRAIN == P) return f_c / elastic_modulus;
    if(ParameterType::CRACKSTRAIN == P) return f_t / elastic_modulus;
    return 0.;
}

unique_ptr<Material> ConcreteExp::get_copy() { return make_unique<ConcreteExp>(*this); }
