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

constexpr unsigned TimberPD::max_iteration = 20u;

const uword TimberPD::sa{0};
const span TimberPD::sb{1, 6};

TimberPD::TimberPD(const unsigned T, vec&& EE, vec&& VV, vec&& SS, vec&& HH, const double R)
    : DataTimberPD{std::forward<vec>(EE), std::forward<vec>(VV), std::forward<vec>(SS), {}, {}, {}, {}, HH(0), HH(1), HH(2), HH(3), HH(4), HH(5), HH(6)}
    , Material3D(T, R) {
    // S(0) = \sigma_{11}^t    S(1) = \sigma_{11}^c
    // S(2) = \sigma_{22}^t    S(3) = \sigma_{22}^c
    // S(4) = \sigma_{33}^t    S(5) = \sigma_{33}^c
    // S(6) = \sigma_{12}^0    S(7) = \sigma_{23}^0    S(8) = \sigma_{13}^0

    proj_a.zeros(6, 6);
    proj_b.zeros(6);

    const auto T1 = 1. / yield_stress(0) / yield_stress(1);
    const auto T2 = 1. / yield_stress(2) / yield_stress(3);
    const auto T3 = 1. / yield_stress(4) / yield_stress(5);

    proj_b(0) = (yield_stress(1) - yield_stress(0)) * (proj_a(0, 0) = T1);
    proj_b(1) = (yield_stress(3) - yield_stress(2)) * (proj_a(1, 1) = T2);
    proj_b(2) = (yield_stress(5) - yield_stress(4)) * (proj_a(2, 2) = T3);

    proj_a(0, 1) = proj_a(1, 0) = -.5 * (T1 + T2 - T3);
    proj_a(1, 2) = proj_a(2, 1) = -.5 * (T2 + T3 - T1);
    proj_a(2, 0) = proj_a(0, 2) = -.5 * (T3 + T1 - T2);
    proj_a(3, 3) = 1. / yield_stress(6) / yield_stress(6);
    proj_a(4, 4) = 1. / yield_stress(7) / yield_stress(7);
    proj_a(5, 5) = 1. / yield_stress(8) / yield_stress(8);
    proj_a *= 2.;

    access::rw(tolerance) = 1E-13;

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
    trial_stiffness = current_stiffness = initial_stiffness = tensor::orthotropic_stiffness(modulus, ratio);

    initial_history.zeros(3);
    initial_history(1) = ini_r_t;
    initial_history(2) = ini_r_c;
    initialize_history(9);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> TimberPD::get_copy() { return make_unique<TimberPD>(*this); }

int TimberPD::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= tolerance) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& k = trial_history(0);
    vec plastic_strain(&trial_history(3), 6, false, true);
    const auto& current_k = current_history(0);

    // plasticity part

    vec principal_stress;    // 3
    mat principal_direction; // 3x3
    if(!eig_sym(principal_stress, principal_direction, tensor::stress::to_tensor(initial_stiffness * (trial_strain - plastic_strain)), "std")) return SUANPAN_FAIL;

    vec principal_tensile_stress(principal_stress.n_elem, fill::zeros), principal_compressive_stress(principal_stress.n_elem, fill::zeros);
    for(auto I = 0llu; I < principal_stress.n_elem; ++I)
        if(principal_stress(I) > 0.) principal_tensile_stress(I) = principal_stress(I);
        else principal_compressive_stress(I) = principal_stress(I);

    mat stiffness_t = transform::eigen_to_tensile_derivative(principal_stress, principal_direction) * initial_stiffness;
    mat stiffness_c = initial_stiffness - stiffness_t;

    const auto trans_mat = transform::compute_jacobian_principal_to_nominal(principal_direction);
    vec sigma_c = trans_mat * principal_compressive_stress;
    const auto trial_sigma_c = sigma_c;

    auto gamma = 0., ref_error = 1.;

    vec incre, residual(7, fill::none);
    mat jacobian(7, 7, fill::none);

    auto counter = 0u;

    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const vec factor_a = proj_a * sigma_c;
        const auto factor_b = dot(proj_b, sigma_c);

        const auto sigma_y = 1. + h * (current_k + gamma);

        const auto f = .5 * dot(factor_a, sigma_c) + factor_b * sigma_y - sigma_y * sigma_y;

        if(1u == counter) {
            if(f < 0.) break;
            gamma = tensor::strain::norm(incre_strain);
        }

        const vec n = factor_a + proj_b * sigma_y;
        const auto norm_n = tensor::strain::norm(n);
        const vec m = n / norm_n;
        const vec unit_m = m % tensor::strain::norm_weight;

        residual(sa) = f;
        residual(sb) = sigma_c + initial_stiffness * gamma * m - trial_sigma_c;

        jacobian(sa, sa) = (factor_b - 2. * sigma_y) * h;
        jacobian(sa, sb) = n.t();

        jacobian(sb, sa) = initial_stiffness * (m + gamma / norm_n * h * (proj_b - m * dot(unit_m, proj_b)));
        jacobian(sb, sb) = eye(6, 6) + gamma / norm_n * initial_stiffness * (proj_a - m * unit_m.t() * proj_a);

        if(!solve(incre, jacobian, residual, solve_opts::refine + solve_opts::equilibrate)) return SUANPAN_FAIL;

        auto error = norm(residual);
        if(1u == counter && error > ref_error) ref_error = error;
        suanpan_debug("Local plasticity iteration error: {:.5E}.\n", error/ref_error);

        if(error <= tolerance * std::max(1., ref_error)) {
            plastic_strain += gamma * m;
            k = current_k + gamma;

            mat left, right(7, 6, fill::zeros);
            right.rows(sb) = stiffness_c;

            if(!solve(left, jacobian, right)) return SUANPAN_FAIL;

            stiffness_c = left.rows(sb);

            break;
        }

        gamma -= incre(sa);
        sigma_c -= incre(sb);
    }

    // damage part

    const vec sigma_t = trans_mat * principal_tensile_stress;

    const auto omega_t = update_damage_t(sigma_t, stiffness_t);
    const auto omega_c = update_damage_c(sigma_c, stiffness_c);

    trial_stress = (1. - omega_t) * sigma_t + (1. - omega_c) * sigma_c;
    trial_stiffness = stiffness_t + stiffness_c;

    return SUANPAN_SUCCESS;
}

double TimberPD::update_damage_t(const vec& sigma_t, mat& stiffness_t) {
    auto& r_t = trial_history(1);

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
    auto& r_c = trial_history(2);

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

int TimberPD::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int TimberPD::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int TimberPD::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

vector<vec> TimberPD::record(const OutputType T) {
    if(T == OutputType::KAPPAP) return {vec{current_history(0)}};
    if(T == OutputType::DT) return {vec{compute_damage_t(current_history(1))}};
    if(T == OutputType::DC) return {vec{compute_damage_c(current_history(2))}};

    return Material3D::record(T);
}

void TimberPD::print() {
    suanpan_info("A 3D Timber Model. doi: 10.1016/j.compstruc.2017.09.010\n");
}
