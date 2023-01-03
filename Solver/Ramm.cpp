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

#include "Ramm.h"
#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/FactoryHelper.hpp>
#include <Solver/Integrator/Integrator.h>

Ramm::Ramm(const unsigned T, const double L, const bool F)
    : Solver(T)
    , arc_length(L)
    , fixed_arc_length(F) {}

int Ramm::analyze() {
    auto& C = get_converger();
    auto& G = get_integrator();
    auto& W = G->get_domain()->get_factory();

    sp_info("-> Current Load Level: {:+.5f}.\n", W->get_trial_load_factor().at(0));

    const auto max_iteration = C->get_max_iteration();

    double t_lambda;

    vec samurai, disp_a, disp_ref;

    // iteration counter
    unsigned counter = 0;

    while(true) {
        // update for nodes and elements
        if(SUANPAN_SUCCESS != G->update_trial_status()) return SUANPAN_FAIL;
        // process modifiers
        if(SUANPAN_SUCCESS != G->process_modifier()) return SUANPAN_FAIL;
        // assemble resistance
        G->assemble_resistance();
        // assemble stiffness
        G->assemble_matrix();
        // process loads
        if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
        // process constraints
        if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;

        // solve ninja
        if(SUANPAN_SUCCESS != G->solve(samurai, G->get_displacement_residual())) return SUANPAN_FAIL;
        // solve reference displacement
        if(SUANPAN_SUCCESS != G->solve(disp_a, G->get_reference_load())) return SUANPAN_FAIL;

        if(const auto n_size = W->get_size(); 0 != W->get_mpc()) {
            mat right, kernel;
            auto& border = W->get_auxiliary_stiffness();
            if(SUANPAN_SUCCESS != G->solve(right, border)) return SUANPAN_FAIL;
            auto& aux_lambda = get_auxiliary_lambda(W);
            if(!solve(aux_lambda, kernel = border.t() * right.head_rows(n_size), border.t() * samurai.head(n_size) - G->get_auxiliary_residual())) return SUANPAN_FAIL;
            samurai -= right * aux_lambda;
            disp_a -= right * solve(kernel, border.t() * disp_a.head_rows(n_size));
        }

        if(0 < counter) t_lambda = -dot(disp_ref, samurai) / dot(disp_ref, disp_a);
        else {
            t_lambda = arc_length / sqrt(dot(disp_a, disp_a) + 1.);

            // check the sign of stiffness for unloading
            if(W->get_stiffness()->sign_det() < 0) t_lambda = -t_lambda;
        }

        // abaqus update
        disp_ref = disp_a;

        samurai += disp_a * t_lambda;

        // avoid machine error accumulation
        G->erase_machine_error(samurai);

        // exit if converged
        if(C->is_converged(counter)) {
            if(!fixed_arc_length) arc_length *= sqrt(max_iteration / static_cast<double>(counter));
            return SUANPAN_SUCCESS;
        }
        // exit if maximum iteration is hit
        if(++counter > max_iteration) {
            if(!fixed_arc_length) arc_length *= .5;
            return SUANPAN_FAIL;
        }

        // update trial displacement
        G->update_from_ninja();
        // update trial load factor
        G->update_trial_load_factor(vec{t_lambda});
        // set time to load factor
        G->update_trial_time(W->get_trial_load_factor().at(0));
        // for tracking
        G->update_load();
        // for tracking
        G->update_constraint();
    }
}

void Ramm::print() { sp_info("A solver using Ramm's arc length method.\n"); }
