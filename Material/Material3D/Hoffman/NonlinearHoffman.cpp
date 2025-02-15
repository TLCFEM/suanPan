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

#include "NonlinearHoffman.h"
#include <Toolbox/tensor.h>

const double NonlinearHoffman::root_two_third = sqrt(2. / 3.);
const uword NonlinearHoffman::sa{0};
const span NonlinearHoffman::sb{1, 6};

NonlinearHoffman::NonlinearHoffman(const unsigned T, vec&& EE, vec&& VV, vec&& SS, const double R)
    : DataNonlinearHoffman{std::move(EE), std::move(VV), std::move(SS)}
    , Material3D(T, R) { transform::hoffman_projection(yield_stress, proj_a, proj_b); }

int NonlinearHoffman::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::orthotropic_stiffness(modulus, ratio);

    elastic_a = initial_stiffness * proj_a;

    initialize_history(7);

    return SUANPAN_SUCCESS;
}

int NonlinearHoffman::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& eqv_strain = trial_history(0);
    const auto& current_eqv_strain = current_history(0);
    vec plastic_strain(&trial_history(1), 6, false, true);

    const vec predictor = (trial_stiffness = initial_stiffness) * (trial_strain - plastic_strain);
    trial_stress = predictor;
    const vec c_stress = .5 * proj_a * initial_stiffness * (current_strain - plastic_strain);

    auto gamma = 0., ref_error = 1.;

    vec7 incre, residual;
    mat77 jacobian;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const vec factor_a = proj_a * trial_stress;
        const vec factor_b = .5 * factor_a + proj_b;
        const vec n_mid = c_stress + factor_b;
        const auto norm_n_mid = root_two_third * tensor::strain::norm(n_mid);
        const auto k = compute_k(eqv_strain = current_eqv_strain + gamma * norm_n_mid);
        const auto f = dot(trial_stress, factor_b) - k * k;

        if(1u == counter && f <= 0.) return SUANPAN_SUCCESS;

        const rowvec dn = two_third / norm_n_mid * (n_mid % tensor::strain::norm_weight).t();
        const auto factor_c = k * compute_dk(eqv_strain);

        residual(sa) = f;
        residual(sb) = trial_stress + gamma * initial_stiffness * n_mid - predictor;

        jacobian(sa, sa) = -2. * factor_c * norm_n_mid;
        jacobian(sa, sb) = factor_a.t() + proj_b.t() - factor_c * gamma * dn * proj_a;
        jacobian(sb, sa) = initial_stiffness * n_mid;
        jacobian(sb, sb) = eye(6, 6) + .5 * gamma * elastic_a;

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local plasticity iteration error: {:.5E}.\n", error);

        if(error < tolerance * ref_error || (inf_norm(residual) < tolerance && counter > 5u)) {
            plastic_strain += gamma * n_mid;

            mat left, right(7, 6, fill::zeros);
            right.rows(sb) = initial_stiffness;

            if(!solve(left, jacobian, right)) return SUANPAN_FAIL;

            trial_stiffness = left.rows(sb);

            return SUANPAN_SUCCESS;
        }

        gamma -= incre(sa);
        trial_stress -= incre(sb);
    }
}

int NonlinearHoffman::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int NonlinearHoffman::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int NonlinearHoffman::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void NonlinearHoffman::print() {
    suanpan_info("A 3D nonlinear hardening model using Hoffman yielding criterion.\n");
}
