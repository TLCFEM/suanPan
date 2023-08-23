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

#include "TimberPD.h"
#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

TimberPD::TimberPD(const unsigned T, vec&& EE, vec&& VV, vec&& SS, vec&& HH, const double R)
    : DataTimberPD{HH(1), HH(2), HH(3), HH(4), HH(5), HH(6)}
    , BilinearHoffman(T, std::forward<vec>(EE), std::forward<vec>(VV), std::forward<vec>(SS), HH(0), R) {
    const auto fill_hill = [&](mat& hill, const double S1, const double S2, const double S3) {
        hill.zeros(6, 6);

        const auto F1 = 2. / S1 / S1;
        const auto F2 = 2. / S2 / S2;
        const auto F3 = 2. / S3 / S3;

        hill(0, 0) = F1;
        hill(1, 1) = F2;
        hill(2, 2) = F3;

        hill(0, 1) = hill(1, 0) = -.5 * (F1 + F2 - F3);
        hill(1, 2) = hill(2, 1) = -.5 * (F2 + F3 - F1);
        hill(2, 0) = hill(0, 2) = -.5 * (F3 + F1 - F2);
        hill(3, 3) = 2. / yield_stress(6) / yield_stress(6);
        hill(4, 4) = 2. / yield_stress(7) / yield_stress(7);
        hill(5, 5) = 2. / yield_stress(8) / yield_stress(8);
    };

    fill_hill(hill_t, yield_stress(0), yield_stress(2), yield_stress(4));
    fill_hill(hill_c, yield_stress(1), yield_stress(3), yield_stress(5));
}

int TimberPD::initialize(const shared_ptr<DomainBase>&) {
    initial_history.resize(9);
    initial_history(7) = ini_r_t;
    initial_history(8) = ini_r_c;
    initialize_history(9);

    return BilinearHoffman::initialize(nullptr);
}

unique_ptr<Material> TimberPD::get_copy() { return make_unique<TimberPD>(*this); }

int TimberPD::update_trial_status(const vec& t_strain) {
    if(SUANPAN_SUCCESS != BilinearHoffman::update_trial_status(t_strain)) return SUANPAN_FAIL;

    if(norm(incre_strain) <= tolerance) return SUANPAN_SUCCESS;

    vec principal_stress;    // 3
    mat principal_direction; // 3x3
    if(!eig_sym(principal_stress, principal_direction, tensor::stress::to_tensor(trial_stress), "std")) return SUANPAN_FAIL;

    mat stiffness_t = transform::eigen_to_tensile_derivative(principal_stress, principal_direction);
    mat stiffness_c = eye(6, 6) - stiffness_t;

    const vec sigma_t = transform::eigen_to_tensile_stress(principal_stress, principal_direction);
    const vec sigma_c = trial_stress - sigma_t;

    const auto omega_t = update_damage_t(sigma_t, stiffness_t);
    const auto omega_c = update_damage_c(sigma_c, stiffness_c);

    trial_stress = (1. - omega_t) * sigma_t + (1. - omega_c) * sigma_c;
    trial_stiffness = (stiffness_t + stiffness_c) * trial_stiffness;

    return SUANPAN_SUCCESS;
}

double TimberPD::update_damage_t(const vec& sigma_t, mat& stiffness_t) {
    auto& r_t = trial_history(7);

    bool new_damage_t = false;
    if(const auto eqv_stress_t = sqrt(.5 * dot(hill_t * sigma_t, sigma_t)); eqv_stress_t > r_t) {
        new_damage_t = true;
        r_t = eqv_stress_t;
    }

    const auto omega_t = compute_damage_t(r_t);
    if(new_damage_t) {
        const auto domega_t = ini_r_t / r_t / r_t * ((b_t * n_t * r_t + n_t) * exp(b_t * (ini_r_t - r_t)) - n_t + 1.);
        const rowvec drdsigma = .5 / r_t * sigma_t.t() * hill_t;
        stiffness_t = ((1. - omega_t) * eye(6, 6) - sigma_t * domega_t * drdsigma) * stiffness_t;
    }
    else stiffness_t *= 1. - omega_t;

    return omega_t;
}

double TimberPD::update_damage_c(const vec& sigma_c, mat& stiffness_c) {
    auto& r_c = trial_history(8);

    bool new_damage_c = false;
    if(const auto eqv_stress_c = sqrt(.5 * dot(hill_c * sigma_c, sigma_c)); eqv_stress_c > r_c) {
        new_damage_c = true;
        r_c = eqv_stress_c;
    }

    const auto omega_c = compute_damage_c(r_c);
    if(new_damage_c) {
        const auto domega_c = m_c * ini_r_c / r_c * omega_c / (r_c - ini_r_c);
        const rowvec drdsigma = .5 / r_c * sigma_c.t() * hill_c;
        stiffness_c = ((1. - omega_c) * eye(6, 6) - sigma_c * domega_c * drdsigma) * stiffness_c;
    }
    else stiffness_c *= 1. - omega_c;

    return omega_c;
}

double TimberPD::compute_damage_t(const double r_t) const { return 1. - ini_r_t / r_t * (1. - n_t + n_t * exp(b_t * (ini_r_t - r_t))); }

double TimberPD::compute_damage_c(const double r_c) const { return b_c * pow(std::max(datum::eps, 1. - ini_r_c / r_c), m_c); }

vector<vec> TimberPD::record(const OutputType T) {
    if(T == OutputType::KAPPAP) return {vec{current_history(0)}};
    if(T == OutputType::DT) return {vec{compute_damage_t(current_history(1))}};
    if(T == OutputType::DC) return {vec{compute_damage_c(current_history(2))}};

    return Material3D::record(T);
}

void TimberPD::print() {
    suanpan_info("A 3D Timber Model. doi: 10.1016/j.compstruc.2017.09.010\n");
}
