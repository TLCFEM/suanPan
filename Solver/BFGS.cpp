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

#include "BFGS.h"

#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>

BFGS::BFGS(const unsigned T, const unsigned MH)
    : Solver(T)
    , max_hist(MH) {}

int BFGS::analyze() {
    auto& C = get_converger();
    auto& G = get_integrator();
    const auto D = C->get_domain().lock();
    auto& W = D->get_factory();

    suanpan_highlight(">> Current Analysis Time: {:.5f}.\n", W->get_trial_time());

    const auto n_size = W->get_size();

    const auto max_iteration = C->get_max_iteration();
    const auto max_storage = 0 == max_hist ? n_size : std::min(max_hist, n_size);

    // iteration counter
    auto counter = 0u;

    vec samurai, residual;

    // clear container
    hist_s.clear();
    hist_y.clear();
    hist_rho.clear();

    wall_clock t_clock;

    while(true) {
        set_step_amplifier(sqrt(max_iteration / (counter + 1.)));

        t_clock.tic();
        // update for nodes and elements
        if(SUANPAN_SUCCESS != G->update_trial_status()) return SUANPAN_FAIL;
        if(SUANPAN_SUCCESS != G->process_modifier()) return SUANPAN_FAIL;
        D->update<Statistics::UpdateStatus>(t_clock.toc());

        t_clock.tic();
        G->assemble_resistance();
        D->update<Statistics::AssembleVector>(t_clock.toc());

        if(0 == counter) {
            t_clock.tic();
            // assemble stiffness for the first iteration
            G->assemble_matrix();
            D->update<Statistics::AssembleMatrix>(t_clock.toc());

            t_clock.tic();
            // process loads and constraints
            if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
            if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;
            D->update<Statistics::ProcessConstraint>(t_clock.toc());

            // indicate the global matrix has been assembled
            G->set_matrix_assembled_switch(true);

            if(0 != W->get_multiplier_size()) {
                suanpan_error("(L-)BFGS solver does not support constraints implemented via the Lagrange multiplier method.\n");
                return SUANPAN_FAIL;
            }

            t_clock.tic();
            // solve the system and commit current displacement increment
            if(SUANPAN_SUCCESS != G->solve(samurai, residual = G->get_force_residual())) return SUANPAN_FAIL;
            D->update<Statistics::SolveSystem>(t_clock.toc());
        }
        else {
            t_clock.tic();
            // process resistance of loads and constraints
            if(SUANPAN_SUCCESS != G->process_load_resistance()) return SUANPAN_FAIL;
            if(SUANPAN_SUCCESS != G->process_constraint_resistance()) return SUANPAN_FAIL;
            D->update<Statistics::ProcessConstraint>(t_clock.toc());

            t_clock.tic();

            // clear temporary factor container
            hist_alpha.clear();
            hist_alpha.reserve(hist_s.size());
            // complete residual increment
            hist_y.back() -= residual = G->get_force_residual();
            // commit current factor after obtaining residual
            hist_rho.emplace_back(dot(hist_s.back(), hist_y.back()));
            // copy current residual to ninja
            samurai = residual;
            // perform two-step recursive loop
            // right side loop
            for(auto J = static_cast<int>(hist_rho.size()) - 1; J >= 0; --J) samurai -= hist_alpha.emplace_back(dot(hist_s[J], samurai)) / hist_rho[J] * hist_y[J];
            // apply the Hessian from the factorization in the first iteration
            samurai = G->solve(samurai);
            // left side loop
            for(size_t I = 0, J = hist_rho.size() - 1; I < hist_rho.size(); ++I, --J) samurai += (hist_alpha[J] - dot(hist_y[I], samurai)) / hist_rho[I] * hist_s[I];

            D->update<Statistics::SolveSystem>(t_clock.toc());
        }

        // commit current displacement increment
        hist_s.emplace_back(samurai);  // complete
        hist_y.emplace_back(residual); // part of residual increment

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
        if(counter > max_storage) {
            hist_s.pop_front();
            hist_y.pop_front();
            hist_rho.pop_front();
        }
    }
}

void BFGS::print() {
    suanpan_info("A (L-)BFGS solver.\n");
}
