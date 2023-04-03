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

#include "BFGS.h"
#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>

BFGS::BFGS(const unsigned T, const unsigned MH)
    : Solver(T)
    , max_hist(std::max(1u, MH)) {}

int BFGS::analyze() {
    auto& C = get_converger();
    auto& G = get_integrator();
    const auto& D = C->get_domain().lock();
    auto& W = D->get_factory();

    suanpan_highlight(">> Current Analysis Time: {:.5f}.\n", W->get_trial_time());

    const auto max_iteration = C->get_max_iteration();

    // iteration counter
    unsigned counter = 0;

    // lambda alias
    auto& aux_lambda = W->modify_auxiliary_lambda();
    vec samurai, residual;

    // clear container
    hist_ninja.clear();
    hist_residual.clear();
    hist_factor.clear();

    auto adjust_for_mpc = [&] {
        if(0 == W->get_mpc()) return SUANPAN_SUCCESS;
        const auto n_size = W->get_size();
        auto& border = W->get_auxiliary_stiffness();
        mat right;
        if(SUANPAN_SUCCESS != G->solve(right, border)) return SUANPAN_FAIL;
        if(!solve(aux_lambda, border.t() * right.head_rows(n_size), border.t() * samurai.head(n_size) - G->get_auxiliary_residual())) return SUANPAN_FAIL;
        samurai -= right * aux_lambda;
        return SUANPAN_SUCCESS;
    };

    while(true) {
        // update for nodes and elements
        if(SUANPAN_SUCCESS != G->update_trial_status()) return SUANPAN_FAIL;
        // process modifiers
        if(SUANPAN_SUCCESS != G->process_modifier()) return SUANPAN_FAIL;
        // assemble resistance
        G->assemble_resistance();

        if(0 == counter) {
            // assemble stiffness for the first iteration
            G->assemble_matrix();
            // process loads and constraints
            if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
            if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;
            // indicate the global matrix has been assembled
            G->set_matrix_assembled_switch(true);
            // solve the system and commit current displacement increment
            if(SUANPAN_SUCCESS != G->solve(samurai, residual = G->get_force_residual())) return SUANPAN_FAIL;
            // deal with mpc
            if(SUANPAN_SUCCESS != adjust_for_mpc()) return SUANPAN_FAIL;
        }
        else {
            // process resistance of loads and constraints
            if(SUANPAN_SUCCESS != G->process_load_resistance()) return SUANPAN_FAIL;
            if(SUANPAN_SUCCESS != G->process_constraint_resistance()) return SUANPAN_FAIL;
            // clear temporary factor container
            alpha.clear();
            alpha.reserve(hist_ninja.size());
            // complete residual increment
            hist_residual.back() -= residual = G->get_force_residual();
            // commit current factor after obtaining residual
            hist_factor.emplace_back(dot(hist_ninja.back(), hist_residual.back()));
            // copy current residual to ninja
            samurai = residual;
            // perform two-step recursive loop
            // right side loop
            for(auto J = static_cast<int>(hist_factor.size()) - 1; J >= 0; --J) {
                // compute and commit alpha
                alpha.emplace_back(dot(hist_ninja[J], samurai) / hist_factor[J]);
                // update ninja
                samurai -= alpha.back() * hist_residual[J];
            }
            // apply the Hessian from the factorization in the first iteration
            samurai = G->solve(samurai);
            // deal with mpc
            if(SUANPAN_SUCCESS != adjust_for_mpc()) return SUANPAN_FAIL;
            // left side loop
            for(size_t I = 0, J = hist_factor.size() - 1; I < hist_factor.size(); ++I, --J) samurai += (alpha[J] - dot(hist_residual[I], samurai) / hist_factor[I]) * hist_ninja[I];
        }

        // commit current displacement increment
        hist_ninja.emplace_back(samurai);     // complete
        hist_residual.emplace_back(residual); // part of residual increment

        // avoid machine error accumulation
        G->erase_machine_error(samurai);

        // exit if converged
        if(C->is_converged(counter)) return G->sync_status(true);
        // exit if maximum iteration is hit
        if(++counter > max_iteration) return SUANPAN_FAIL;

        // update internal variable
        G->update_internal(samurai);
        // update trial status for factory
        G->update_from_ninja();
        // for tracking
        G->update_load();
        // for tracking multiplier
        G->update_constraint();

        // fast handling for linear elastic case
        // sync status using newly computed increment across elements and nodes
        // this may just call predictor or call corrector
        if(D->get_attribute(ModalAttribute::LinearSystem)) return G->sync_status(false);

        // check if the maximum record number is hit (L-BFGS)
        if(counter > max_hist) {
            hist_ninja.pop_front();
            hist_residual.pop_front();
            hist_factor.pop_front();
        }
    }
}

void BFGS::print() {
    suanpan_info("A (L-)BFGS solver.\n");
}
