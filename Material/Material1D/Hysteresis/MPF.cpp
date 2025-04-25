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

#include "MPF.h"
#include <Toolbox/utility.h>

MPF::MPF(const unsigned T, const double E, const double Y, const double H, const double R, const double B1, const double B2, const double B3, const double B4, const bool ISO, const bool CON, const double D)
    : DataMPF{fabs(E), H, fabs(Y), fabs(R), fabs(B1), fabs(B2), fabs(B3), fabs(B4), ISO, CON}
    , Material1D(T, D) {}

int MPF::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(7);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> MPF::get_copy() { return make_unique<MPF>(*this); }

int MPF::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& reverse_stress = trial_history(0);
    auto& reverse_strain = trial_history(1);
    auto& inter_stress = trial_history(2);
    auto& inter_strain = trial_history(3);
    auto& pre_inter_strain = trial_history(4);
    auto& max_strain = trial_history(5);
    auto& load_sign = trial_history(6);

    auto shift_stress = 0.;
    if(isotropic_hardening) {
        shift_stress = std::max(0., A3 * yield_stress * (max_strain / yield_strain - A4));
        max_strain = std::max(max_strain, fabs(trial_strain(0)));
    }

    if(const auto trial_load_sign = suanpan::sign(incre_strain(0)); !suanpan::approx_equal(trial_load_sign, load_sign)) {
        if(!suanpan::approx_equal(load_sign, 0.)) {
            reverse_stress = current_stress(0);
            reverse_strain = current_strain(0);
            pre_inter_strain = inter_strain;
            inter_strain = yield_stress * hardening_ratio - yield_stress - shift_stress;
            if(trial_load_sign > 0.) inter_strain = -inter_strain;
            inter_strain = (inter_strain + elastic_modulus * reverse_strain - reverse_stress) / (elastic_modulus - hardening_ratio * elastic_modulus);
            inter_stress = elastic_modulus * (inter_strain - reverse_strain) + reverse_stress;
        }
        else if(trial_load_sign > 0.) {
            inter_stress = yield_stress;
            inter_strain = yield_strain;
        }
        else {
            inter_stress = -yield_stress;
            inter_strain = -yield_strain;
        }
        load_sign = trial_load_sign;
    }

    auto radius = R0;
    if(!constant_radius) {
        // update radius
        const auto xi = fabs(reverse_strain - pre_inter_strain) / yield_strain;
        radius -= A1 * xi / (A2 + xi);
    }

    const auto gap_strain = inter_strain - reverse_strain;
    const auto gap_stress = inter_stress - reverse_stress;
    const auto normal_strain = std::max(datum::eps, (trial_strain(0) - reverse_strain) / gap_strain);
    const auto factor_a = 1. + pow(normal_strain, radius);
    const auto factor_b = (1. - hardening_ratio) * pow(factor_a, -1. / radius);

    trial_stress = (hardening_ratio + factor_b) * normal_strain * gap_stress + reverse_stress;
    trial_stiffness = gap_stress / gap_strain * (hardening_ratio + factor_b / factor_a);

    suanpan_assert([&] { if(!trial_stress.is_finite() || !trial_stiffness.is_finite()) throw std::invalid_argument("infinite number detected"); });

    return SUANPAN_SUCCESS;
}

int MPF::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int MPF::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    current_history = trial_history;
    return SUANPAN_SUCCESS;
}

int MPF::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    trial_history = current_history;
    return SUANPAN_SUCCESS;
}

void MPF::print() {
    suanpan_info("A Menegotto-Pinto-Filippou model with initial stiffness {:.3E} and yield stress {:.3E} with isotropic hardening {} and Bauschinger effect {}.\n", elastic_modulus, yield_stress, isotropic_hardening ? "enabled" : "disabled", constant_radius ? "disabled" : "enabled");
    Material1D::print();
}
