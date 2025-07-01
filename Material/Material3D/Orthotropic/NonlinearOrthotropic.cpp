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
const span NonlinearOrthotropic::sb{1, 6};

NonlinearOrthotropic::NonlinearOrthotropic(const unsigned T, const OrthotropicType TP, vec&& EE, vec&& VV, vec&& SS, const double R)
    : DataNonlinearOrthotropic{std::move(EE), std::move(VV), std::move(SS)}
    , Material3D(T, R) {
    switch(TP) {
    case OrthotropicType::Hoffman:
        transform::hoffman_projection(yield_stress, proj_p, proj_q);
        break;
    case OrthotropicType::TsaiWu:
        transform::tsai_wu_projection(yield_stress, proj_p, proj_q);
        break;
    }

    // a quite arbitrary and empirical value that does not necessarily imply any physical meaning
    // just to provide a threshold to switch between two integration methods
    reference_strain = tensor::strain::norm(vec{.5 * (yield_stress(0) + yield_stress(1)), .5 * (yield_stress(2) + yield_stress(3)), .5 * (yield_stress(4) + yield_stress(5)), yield_stress(6), yield_stress(7), yield_stress(8)} / modulus);
}

int NonlinearOrthotropic::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::orthotropic_stiffness(modulus, ratio);

    elastic_p = initial_stiffness * proj_p;

    initialize_history(7);

    return SUANPAN_SUCCESS;
}

int NonlinearOrthotropic::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    const auto norm_incre_strain = tensor::strain::norm(incre_strain);

    if(norm_incre_strain <= datum::eps) return SUANPAN_SUCCESS;

    return norm_incre_strain < reference_strain ? trapezoidal_return() : euler_return();
}

int NonlinearOrthotropic::trapezoidal_return() {
    trial_history = current_history;
    const auto& current_ep = current_history(0);
    auto& ep = trial_history(0);
    vec plastic_strain(&trial_history(1), 6, false, true);

    const vec6 base_s = (trial_stiffness = initial_stiffness) * (current_strain - plastic_strain);
    const vec6 incre_s = initial_stiffness * incre_strain;
    const vec6 predictor_s = trial_stress = base_s + incre_s;

    const auto current_k = compute_k(current_ep);

    // elastic region
    if(ortho_inner(predictor_s) <= current_k * current_k) return SUANPAN_SUCCESS;

    vec6 onset_s, onset_n;
    mat66 dnds;

    const auto inside_surface = ortho_inner(base_s) < current_k * current_k;
    if(inside_surface) {
        // elastic loading to plastic
        // need to find the onset point
        const auto onset_r = ridders([&](const double r) { return ortho_inner(onset_s = base_s + r * incre_s) - current_k * current_k; }, 0., 1., tolerance);
        onset_n = proj_p * onset_s + proj_q;
        dnds = proj_p * (onset_r * eye(6, 6) - onset_r / dot(onset_n, incre_s) * incre_s * onset_n.t());
    }
    else onset_n = proj_p * (onset_s = base_s) + proj_q;

    auto gamma{0.};

    vec7 incre(fill::none), residual(fill::none);
    mat77 jacobian(fill::none);

    auto counter{0u};
    auto try_bisection{false};
    while(true) {
        if(max_iteration == ++counter) {
            try_bisection = true;

            const vec6 const_a = onset_n + proj_q;
            const vec6 const_b = initial_stiffness * const_a;

            const auto approx_update = [&](const double gm) {
                gamma = gm;
                const vec6 approx_a = proj_p * (trial_stress = solve(eye(6, 6) + .5 * gamma * elastic_p, predictor_s - .5 * gamma * const_b));
                const auto approx_k = compute_k(ep = current_ep + .5 * gamma * root_two_third * tensor::strain::norm(approx_a + const_a));
                return .5 * dot(trial_stress, approx_a) + dot(trial_stress, proj_q) - approx_k * approx_k;
            };

            ridders_guess(approx_update, 0., .25 * tensor::strain::norm(incre_strain) / tensor::strain::norm(proj_p * predictor_s + proj_q), tolerance);
        }

        const vec6 n = proj_p * trial_stress + proj_q;
        const vec6 n_mid = .5 * (n + onset_n);
        const auto norm_n_mid = root_two_third * tensor::strain::norm(n_mid);
        const rowvec6 pnn = .5 * two_third / norm_n_mid * (n_mid % tensor::strain::norm_weight).t(); // part of some derivative
        const auto k = compute_k(ep = current_ep + gamma * norm_n_mid);

        residual(sa) = ortho_inner(trial_stress) - k * k;
        residual(sb) = trial_stress + gamma * initial_stiffness * n_mid - predictor_s;

        const auto pk = -2. * k * compute_dk(ep);

        jacobian(sa, sa) = pk * norm_n_mid;
        jacobian(sa, sb) = n.t() + pk * gamma * pnn * proj_p;
        jacobian(sb, sa) = initial_stiffness * n_mid;
        jacobian(sb, sb) = eye(6, 6) + .5 * gamma * elastic_p;

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        suanpan_debug("Local plasticity iteration (2nd) error: {:.5E}.\n", error);

        if(error < tolerance * (1. + yield_stress.max()) || (inf_norm(residual) < tolerance && counter > 5u) || try_bisection) {
            plastic_strain += gamma * n_mid;

            mat::fixed<7, 6> left(fill::none), right(fill::zeros);
            if(inside_surface) {
                right.row(sa) = -pk * gamma * pnn * dnds;
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

int NonlinearOrthotropic::euler_return() {
    trial_history = current_history;
    const auto& current_ep = current_history(0);
    auto& ep = trial_history(0);
    vec plastic_strain(&trial_history(1), 6, false, true);

    const vec6 predictor_s = trial_stress = (trial_stiffness = initial_stiffness) * (trial_strain - plastic_strain);

    // elastic region
    if(ortho_inner(predictor_s) <= std::pow(compute_k(current_ep), 2.)) return SUANPAN_SUCCESS;

    auto gamma{0.};

    vec7 incre(fill::none), residual(fill::none);
    mat77 jacobian(fill::none);

    auto counter{0u};
    auto try_bisection{false};
    while(true) {
        if(max_iteration == ++counter) {
            try_bisection = true;

            const vec6 elastic_q = initial_stiffness * proj_q;

            const auto approx_update = [&](const double gm) {
                gamma = gm;
                const vec6 approx_a = proj_p * (trial_stress = solve(eye(6, 6) + gamma * elastic_p, predictor_s - gamma * elastic_q));
                const auto approx_k = compute_k(ep = current_ep + gamma * root_two_third * tensor::strain::norm(approx_a + proj_q));
                return .5 * dot(trial_stress, approx_a) + dot(trial_stress, proj_q) - approx_k * approx_k;
            };

            ridders_guess(approx_update, 0., .25 * tensor::strain::norm(incre_strain) / tensor::strain::norm(proj_p * predictor_s + proj_q), tolerance);
        }

        const vec6 n = proj_p * trial_stress + proj_q;
        const auto norm_n = root_two_third * tensor::strain::norm(n);
        const rowvec6 pnn = two_third / norm_n * (n % tensor::strain::norm_weight).t(); // part of some derivative
        const auto k = compute_k(ep = current_ep + gamma * norm_n);

        residual(sa) = ortho_inner(trial_stress) - k * k;
        residual(sb) = trial_stress + gamma * initial_stiffness * n - predictor_s;

        const auto pk = -2. * k * compute_dk(ep);

        jacobian(sa, sa) = pk * norm_n;
        jacobian(sa, sb) = n.t() + pk * gamma * pnn * proj_p;
        jacobian(sb, sa) = initial_stiffness * n;
        jacobian(sb, sb) = eye(6, 6) + gamma * elastic_p;

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        suanpan_debug("Local plasticity iteration (1st) error: {:.5E}.\n", error);

        if(error < tolerance * (1. + yield_stress.max()) || (inf_norm(residual) < tolerance && counter > 5u) || try_bisection) {
            plastic_strain += gamma * n;

            mat::fixed<7, 6> left(fill::none), right(fill::zeros);
            right.rows(sb) = initial_stiffness;

            if(!solve(left, jacobian, right, solve_opts::equilibrate)) return SUANPAN_FAIL;

            trial_stiffness = left.rows(sb);

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

void NonlinearOrthotropic::print() { suanpan_info("A 3D nonlinear hardening model using an orthotropic yielding criterion.\n"); }

NonlinearHill::NonlinearHill(const unsigned T, vec&& EE, vec&& VV, vec&& S, const double R)
    : NonlinearOrthotropic(T, OrthotropicType::Hoffman, std::move(EE), std::move(VV), vec{S(0), S(0), S(1), S(1), S(2), S(2), S(3), S(4), S(5)}, R) {}
