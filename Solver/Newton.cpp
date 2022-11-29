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

#include "Newton.h"
#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/FactoryHelper.hpp>
#include <Solver/Integrator/Integrator.h>

Newton::Newton(const unsigned T, const bool IS)
    : Solver(T)
    , initial_stiffness(IS) {}

int Newton::analyze() {
    auto& C = get_converger();
    auto& G = get_integrator();
    auto& W = G->get_domain().lock()->get_factory();

    suanpan_info("current analysis time: %.5f.\n", W->get_trial_time());

    const auto max_iteration = C->get_max_iteration();

    // iteration counter
    unsigned counter = 0;

    // ninja alias
    auto& ninja = get_ninja(W);

    vec pre_ninja;

    auto aitken = false;

    while(true) {
        // update for nodes and elements
        if(SUANPAN_SUCCESS != G->update_trial_status()) return SUANPAN_FAIL;
        // process modifiers
        if(SUANPAN_SUCCESS != G->process_modifier()) return SUANPAN_FAIL;
        // assemble resistance
        G->assemble_resistance();

        if(initial_stiffness && counter != 0) {
            // some loads may have resistance
            if(SUANPAN_SUCCESS != G->process_load_resistance()) return SUANPAN_FAIL;
            // some constraints may have resistance
            if(SUANPAN_SUCCESS != G->process_constraint_resistance()) return SUANPAN_FAIL;
        }
        else {
            // first iteration
            // assemble stiffness
            G->assemble_matrix();
            // process loads
            if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
            // process constraints
            if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;
        }

        // call solver
        auto flag = G->solve(ninja, G->get_force_residual());

        suanpan_debug([&] {
            if(!ninja.is_finite()) {
                suanpan_fatal("infinite number detected.\n");
                flag = SUANPAN_FAIL;
            }
        });

        // make sure system solver succeeds
        if(SUANPAN_SUCCESS != flag) {
            suanpan_error("system solver returns error code %d.\n", flag);
            return flag;
        }

        // deal with mpc
        if(0 != W->get_mpc()) {
            const auto n_size = W->get_size();
            auto& border = W->get_auxiliary_stiffness();
            mat right;
            if(SUANPAN_SUCCESS != G->solve(right, border)) return SUANPAN_FAIL;
            auto& aux_lambda = get_auxiliary_lambda(W);
            if(!solve(aux_lambda, border.t() * right.head_rows(n_size), border.t() * ninja.head_rows(n_size) - G->get_auxiliary_residual())) return SUANPAN_FAIL;
            ninja -= right * aux_lambda;
        }

        if(initial_stiffness) {
            if(!aitken) {
                aitken = true;
                pre_ninja = ninja;
            }
            else {
                aitken = false;
                const vec diff_ninja = pre_ninja - ninja;
                ninja *= dot(pre_ninja, diff_ninja) / dot(diff_ninja, diff_ninja);
            }
        }

        // avoid machine error accumulation
        G->erase_machine_error();

        // exit if converged
        if(C->is_converged(counter)) return SUANPAN_SUCCESS;
        // exit if maximum iteration is hit
        if(++counter > max_iteration) return SUANPAN_FAIL;

        // update internal variable
        G->update_internal(ninja);
        // update trial status for factory
        G->update_from_ninja(ninja);
        // for tracking
        G->update_load();
        // for tracking
        G->update_constraint();
    }
}

void Newton::print() { suanpan_info("A solver based on Newton--Raphson iteration method %s", initial_stiffness ? "using initial stiffness for each substep.\n" : ".\n"); }
