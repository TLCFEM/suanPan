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

#include "NonlocalBilinearJ2.h"

#include <Toolbox/tensor.h>

constexpr double NonlocalBilinearJ2::two_third = 2. / 3.;
const double NonlocalBilinearJ2::root_two_third = std::sqrt(two_third);
const mat NonlocalBilinearJ2::unit_dev_tensor = tensor::unit_deviatoric_tensor4();
const uvec NonlocalBilinearJ2::UD{0, 1, 2, 3, 4, 5};
const uvec NonlocalBilinearJ2::DD{6};

NonlocalBilinearJ2::NonlocalBilinearJ2(const unsigned T, const double E, const double V, const double Y, const double PE, const double ER, const double H, const double B, const double R)
    : DataNonlocalBilinearJ2{.elastic_modulus = std::fabs(E), .poissons_ratio = V, .yield_stress = std::fabs(Y), .plastic_strain_threshold = PE, .evolution_rate = std::fabs(ER), .hardening_ratio = H, .beta = std::clamp(B, 0., 1.)}
    , Material3D(T, R) {}

int NonlocalBilinearJ2::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio, nonlocal_size());

    initialize_history(13);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> NonlocalBilinearJ2::unique_copy() { return std::make_unique<NonlocalBilinearJ2>(*this); }

double NonlocalBilinearJ2::get(const Parameter P) const { return MaterialProperty(elastic_modulus, poissons_ratio)(P); }

int NonlocalBilinearJ2::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& eqv_plastic_strain = trial_history(0);
    vec back_stress(&trial_history(1), 6, false, true);
    vec plastic_strain(&trial_history(7), 6, false, true);

    trial_stiffness = initial_stiffness;
    trial_stress.head(6) = trial_stiffness(UD, UD) * (trial_strain.head(6) - plastic_strain);

    const vec shifted_stress = tensor::dev(trial_stress.head(6)) - back_stress;

    const auto norm_shifted_stress = tensor::stress::norm(shifted_stress);

    const auto yield_surf = yield_stress + isotropic_modulus * eqv_plastic_strain;
    const auto yield_func = norm_shifted_stress - root_two_third * std::max(yield_surf, 0.);
    const auto yield_flag = yield_func > 0.;

    rowvec ppepe;

    if(yield_flag) {
        const auto tmp_a = double_shear + two_third * (kinematic_modulus + (yield_surf > 0. ? isotropic_modulus : 0.));
        const auto gamma = yield_func / tmp_a;

        auto tmp_b = double_shear * gamma / norm_shifted_stress;
        trial_stress.head(6) -= tmp_b * shifted_stress;

        const vec unit_n = shifted_stress / norm_shifted_stress;

        eqv_plastic_strain += root_two_third * gamma;
        back_stress += two_third * kinematic_modulus * gamma * unit_n;
        plastic_strain += gamma * unit_n % tensor::stress::norm_weight;

        ppepe = root_two_third * double_shear / tmp_a * unit_n.t();

        tmp_b *= double_shear;
        trial_stiffness(UD, UD) += (tmp_b - square_double_shear / tmp_a) * unit_n * unit_n.t() - tmp_b * unit_dev_tensor;
    }

    if(eqv_plastic_strain > plastic_strain_threshold) {
        trial_stress(6) = 1. - std::exp(evolution_rate * (plastic_strain_threshold - eqv_plastic_strain));
        if(yield_flag) trial_stiffness(DD, UD) = (1. - trial_stress(6)) * evolution_rate * ppepe;
    }

    trial_stiffness(UD, DD) = -trial_stress(UD);
    trial_stress(UD) *= 1. - trial_strain(6);
    trial_stiffness(UD, UD) *= 1. - trial_strain(6);

    return SUANPAN_SUCCESS;
}

void NonlocalBilinearJ2::print() {
    suanpan_info("A 3D bilinear hardening model.\n");
}
