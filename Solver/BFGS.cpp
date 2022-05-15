/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
#include <Domain/FactoryHelper.hpp>
#include <Solver/Integrator/Integrator.h>

BFGS::BFGS(const unsigned T, const unsigned MH)
    : Solver(T)
    , max_hist(std::max(1u, MH)) {}

int BFGS::analyze() {
    auto& C = get_converger();
    auto& G = get_integrator();
    const auto& D = C->get_domain().lock();
    auto& W = D->get_factory();

    const auto max_iteration = C->get_max_iteration();

    suanpan_info("current analysis time: %.5f.\n", W->get_trial_time());

    // iteration counter
    unsigned counter = 0;

    // ninja alias
    auto& ninja = get_ninja(W);
    // lambda alias
    auto& aux_lambda = get_auxiliary_lambda(W);

    // clear container
    hist_ninja.clear();
    hist_residual.clear();
    hist_factor.clear();

    while(true) {
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
            // commit current residual
            hist_residual.emplace_back(G->get_force_residual());
            // solve the system and commit current displacement increment
            if(SUANPAN_SUCCESS != G->solve(ninja, hist_residual.back())) return SUANPAN_FAIL;
            // deal with mpc
            if(0 != W->get_mpc()) {
                const auto n_size = W->get_size();
                auto& border = W->get_auxiliary_stiffness();
                mat right;
                if(SUANPAN_SUCCESS != G->solve(right, border)) return SUANPAN_FAIL;
                if(!solve(aux_lambda, border.t() * right.head_rows(n_size), border.t() * ninja.head_rows(n_size) - G->get_auxiliary_residual())) return SUANPAN_FAIL;
                ninja -= right * aux_lambda;
            }
        }
        else {
            // process resistance of loads and constraints
            if(SUANPAN_SUCCESS != G->process_load_resistance()) return SUANPAN_FAIL;
            if(SUANPAN_SUCCESS != G->process_constraint_resistance()) return SUANPAN_FAIL;
            // clear temporary factor container
            alpha.clear();
            // commit current residual
            hist_residual.emplace_back(G->get_force_residual());
            // copy current residual to ninja
            ninja = hist_residual.back();
            // perform two-step recursive loop
            // right side loop
            for(size_t I = 0, J = hist_factor.size() - 1; I < hist_factor.size(); ++I, --J) {
                // compute and commit alpha
                alpha.emplace_back(dot(hist_ninja[J], ninja) / hist_factor[J]);
                // update ninja
                ninja -= alpha.back() * hist_residual[J];
            }
            // apply the Hessian from the factorization in the first iteration
            ninja = G->solve(ninja);
            // deal with mpc
            if(0 != W->get_mpc()) {
                const auto n_size = W->get_size();
                auto& border = W->get_auxiliary_stiffness();
                mat right;
                if(SUANPAN_SUCCESS != G->solve(right, border)) return SUANPAN_FAIL;
                if(!solve(aux_lambda, border.t() * right.head_rows(n_size), border.t() * ninja.head_rows(n_size) - G->get_auxiliary_residual())) return SUANPAN_FAIL;
                ninja -= right * aux_lambda;
            }
            // left side loop
            for(size_t I = 0, J = hist_factor.size() - 1; I < hist_factor.size(); ++I, --J) ninja += (alpha[J] - dot(hist_residual[I], ninja) / hist_factor[I]) * hist_ninja[I];
        }
        // commit current displacement increment
        hist_ninja.emplace_back(ninja);
        // commit current factor after obtaining ninja and residual
        hist_factor.emplace_back(dot(hist_ninja.back(), hist_residual.back()));

        // avoid machine error accumulation
        G->erase_machine_error();
        // update internal variable
        G->update_internal(ninja);
        // update trial status for factory
        G->update_trial_displacement(ninja);
        // for tracking
        G->update_load();
        // for tracking multiplier
        G->update_constraint();
        // update for nodes and elements
        if(SUANPAN_SUCCESS != G->update_trial_status()) return SUANPAN_FAIL;

        // exit if converged
        if(C->is_converged()) return SUANPAN_SUCCESS;
        // exit if maximum iteration is hit
        if(++counter > max_iteration) return SUANPAN_FAIL;

        // check if the maximum record number is hit (L-BFGS)
        if(counter > max_hist) {
            hist_ninja.pop_front();
            hist_residual.pop_front();
            hist_factor.pop_front();
        }
    }
}

void BFGS::print() { suanpan_info("(L-)BFGS.\n"); }
