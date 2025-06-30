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

#include "NonlinearOrthotropic.h"

#include <Toolbox/ridders.hpp>
#include <Toolbox/tensor.h>

const double NonlinearOrthotropic::root_two_third = std::sqrt(2. / 3.);
const uword NonlinearOrthotropic::sa{0};
const span NonlinearOrthotropic::sb{1, 6};

NonlinearOrthotropic::NonlinearOrthotropic(const unsigned T, const OrthotropicType TP, vec&& EE, vec&& VV, vec&& SS, const double R)
    : DataNonlinearOrthotropic{std::move(EE), std::move(VV), std::move(SS)}
    , Material3D(T, R) {
    switch(TP) {
    case OrthotropicType::Hoffman:
        transform::hoffman_projection(yield_stress, proj_a, proj_b);
        break;
    case OrthotropicType::TsaiWu:
        transform::tsai_wu_projection(yield_stress, proj_a, proj_b);
        break;
    }
}

int NonlinearOrthotropic::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::orthotropic_stiffness(modulus, ratio);

    elastic_a = initial_stiffness * proj_a;

    initialize_history(7);

    return SUANPAN_SUCCESS;
}

int NonlinearOrthotropic::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& eqv_strain = trial_history(0);
    const auto& current_eqv_strain = current_history(0);
    vec plastic_strain(&trial_history(1), 6, false, true);

    const vec6 base_stress = (trial_stiffness = initial_stiffness) * (current_strain - plastic_strain);
    const vec6 offset_stress = initial_stiffness * incre_strain;
    const vec6 predictor = base_stress + offset_stress;
    trial_stress = predictor;

    const auto current_k = compute_k(current_eqv_strain);

    if(.5 * dot(trial_stress, proj_a * trial_stress) + dot(trial_stress, proj_b) <= current_k * current_k) return SUANPAN_SUCCESS;

    vec6 onset_stress, onset_n;
    mat66 dnds;

    const auto inside_surface = .5 * dot(base_stress, proj_a * base_stress) + dot(base_stress, proj_b) < current_k * current_k;
    if(inside_surface) {
        // elastic loading to plastic
        const auto compute_onset = [&](const double r) {
            onset_stress = base_stress + r * offset_stress;
            return .5 * dot(onset_stress, proj_a * onset_stress) + dot(onset_stress, proj_b) - current_k * current_k;
        };

        const auto onset_r = ridders(compute_onset, 0., 1., tolerance);
        onset_n = proj_a * onset_stress + proj_b;
        dnds = proj_a * (onset_r * eye(6, 6) - onset_r / dot(onset_n, offset_stress) * offset_stress * onset_n.t());
    }
    else onset_n = proj_a * (onset_stress = base_stress) + proj_b;

    auto gamma = 0.;

    vec7 incre(fill::none), residual(fill::none);
    mat77 jacobian(fill::none);

    auto counter = 0u;
    auto try_bisection = false;
    while(true) {
        if(max_iteration == ++counter) {
            try_bisection = true;

            const auto approx_update = [&](const double gm) {
                gamma = gm;
                const vec6 approx_a = proj_a * (trial_stress = solve(eye(6, 6) + .5 * gamma * elastic_a, predictor - .5 * gamma * initial_stiffness * (onset_n + proj_b)));
                const auto approx_k = compute_k(eqv_strain = current_eqv_strain + .5 * gamma * root_two_third * tensor::strain::norm(approx_a + proj_b + onset_n));
                return .5 * dot(trial_stress, approx_a) + dot(trial_stress, proj_b) - approx_k * approx_k;
            };

            ridders_guess(approx_update, 0., .25 * tensor::strain::norm(incre_strain) / tensor::strain::norm(proj_a * predictor + proj_b), tolerance);
        }

        const vec6 factor_a = proj_a * trial_stress;
        const vec6 n = factor_a + proj_b;
        const vec6 n_mid = .5 * (n + onset_n);
        const auto norm_n_mid = root_two_third * tensor::strain::norm(n_mid);
        const auto k = compute_k(eqv_strain = current_eqv_strain + gamma * norm_n_mid);
        const rowvec6 pnn = .5 * two_third / norm_n_mid * (n_mid % tensor::strain::norm_weight).t(); // part of some derivative

        residual(sa) = .5 * dot(trial_stress, factor_a) + dot(trial_stress, proj_b) - k * k;
        residual(sb) = trial_stress + gamma * initial_stiffness * n_mid - predictor;

        const auto factor_b = -2. * k * compute_dk(eqv_strain);

        jacobian(sa, sa) = factor_b * norm_n_mid;
        jacobian(sa, sb) = n.t() + factor_b * gamma * pnn * proj_a;
        jacobian(sb, sa) = initial_stiffness * n_mid;
        jacobian(sb, sb) = eye(6, 6) + .5 * gamma * elastic_a;

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        suanpan_debug("Local plasticity iteration error: {:.5E}.\n", error);

        if(error < tolerance * (1. + yield_stress.max()) || (inf_norm(residual) < tolerance && counter > 5u) || try_bisection) {
            plastic_strain += gamma * n_mid;

            mat::fixed<7, 6> left(fill::none), right(fill::zeros);
            if(inside_surface) {
                right.row(sa) = -factor_b * gamma * pnn * dnds;
                right.rows(sb) = eye(6, 6) - .5 * gamma * initial_stiffness * dnds;
            }
            else right.rows(sb) = eye(6, 6);

            if(!solve(left, jacobian, right, solve_opts::equilibrate)) return SUANPAN_FAIL;

            trial_stiffness = left.rows(sb) * initial_stiffness;

            return SUANPAN_SUCCESS;
        }

        gamma = std::max(0., gamma - incre(sa));
        trial_stress -= incre(sb);
    }
}

int NonlinearOrthotropic::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int NonlinearOrthotropic::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int NonlinearOrthotropic::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void NonlinearOrthotropic::print() { suanpan_info("A 3D nonlinear hardening model using orthotropic yielding criterion.\n"); }

NonlinearHill::NonlinearHill(const unsigned T, vec&& EE, vec&& VV, vec&& S, const double R)
    : NonlinearOrthotropic(T, OrthotropicType::Hoffman, std::move(EE), std::move(VV), vec{S(0), S(0), S(1), S(1), S(2), S(2), S(3), S(4), S(5)}, R) {}
