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

#include "ConcreteTsai.h"
#include <Toolbox/utility.h>

pod2 ConcreteTsai::compute_compression_initial_reverse() const {
    pod2 response;

    response[1] = compute_compression_backbone(response[0] = middle_point * c_strain)[0];

    return response;
}

pod2 ConcreteTsai::compute_tension_initial_reverse() const {
    pod2 response;

    response[1] = compute_tension_backbone(response[0] = middle_point * t_strain)[0];

    return response;
}

pod2 ConcreteTsai::compute_compression_backbone(const double n_strain) const {
    pod2 response;

    const auto normal_strain = std::max(datum::eps, n_strain / c_strain);
    const auto tmp_a = pow(normal_strain, c_n);
    const auto tmp_b = c_n == 1. ? 1. + (c_m - 1. + log(normal_strain)) * normal_strain : 1. + (c_m - c_n / (c_n - 1.)) * normal_strain + tmp_a / (c_n - 1.);
    response[0] = c_stress * c_m * normal_strain / tmp_b;
    response[1] = initial_stiffness(0) * (1. - tmp_a) / tmp_b / tmp_b;

    return response;
}

pod2 ConcreteTsai::compute_tension_backbone(const double n_strain) const {
    pod2 response;

    const auto normal_strain = std::max(datum::eps, n_strain / t_strain);
    const auto tmp_a = pow(normal_strain, t_n);
    const auto tmp_b = t_n == 1. ? 1. + (t_m - 1. + log(normal_strain)) * normal_strain : 1. + (t_m - t_n / (t_n - 1.)) * normal_strain + tmp_a / (t_n - 1.);
    response[0] = t_stress * t_m * normal_strain / tmp_b;
    response[1] = initial_stiffness(0) * (1. - tmp_a) / tmp_b / tmp_b;

    return response;
}

double ConcreteTsai::compute_compression_residual(const double reverse_c_strain, const double reverse_c_stress) const { return std::min(0., reverse_c_strain - reverse_c_stress * (reverse_c_strain / c_strain + .57) / (reverse_c_stress / c_strain + .57 * initial_stiffness(0))); }

double ConcreteTsai::compute_tension_residual(const double reverse_t_strain, const double reverse_t_stress) const { return std::max(0., reverse_t_strain - reverse_t_stress * (reverse_t_strain / t_strain + .67) / (reverse_t_stress / t_strain + .67 * initial_stiffness(0))); }

ConcreteTsai::ConcreteTsai(const unsigned T, const double E, const double CS, const double TS, const double NC, const double NT, const double MP, const double CE, const double TE, const double R)
    : SimpleHysteresis(T, MP, R)
    , elastic_modulus(fabs(E))
    , c_stress(-perturb(fabs(CS)))
    , c_strain(-perturb(fabs(CE)))
    , c_m(elastic_modulus * c_strain / c_stress)
    , c_n(std::max(perturb(1.), NC))
    , t_stress(perturb(fabs(TS)))
    , t_strain(perturb(fabs(TE)))
    , t_m(elastic_modulus * t_strain / t_stress)
    , t_n(std::max(perturb(1.), NT)) {}

int ConcreteTsai::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    return SUANPAN_SUCCESS;
}

double ConcreteTsai::get_parameter(const ParameterType P) const {
    if(ParameterType::ELASTICMODULUS == P) return initial_stiffness(0);
    if(ParameterType::PEAKSTRAIN == P) return c_strain;
    if(ParameterType::CRACKSTRAIN == P) return t_strain;
    return 0.;
}

unique_ptr<Material> ConcreteTsai::get_copy() { return make_unique<ConcreteTsai>(*this); }
