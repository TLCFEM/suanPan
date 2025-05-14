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

#include "Newton.h"

#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>

Newton::Newton(const unsigned T, const bool IS)
    : Solver(T)
    , initial_stiffness(IS) {}

int Newton::analyze() {
    auto& C = get_converger();
    auto& G = get_integrator();
    const auto D = G->get_domain();
    auto& W = D->get_factory();

    suanpan_highlight(">> Current Analysis Time: {:.5f}.\n", W->get_trial_time());

    const auto max_iteration = C->get_max_iteration();

    // iteration counter
    auto counter = 0u;

    inner_product = datum::eps;

    vec samurai, pre_samurai;

    auto aitken = false;

    wall_clock t_clock;

    while(true) {
        set_step_amplifier(sqrt(max_iteration / (counter + 1.)));

        t_clock.tic();
        if(SUANPAN_SUCCESS != G->update_trial_status()) return SUANPAN_FAIL;
        if(SUANPAN_SUCCESS != G->process_modifier()) return SUANPAN_FAIL;
        D->update<Statistics::UpdateStatus>(t_clock.toc());

        t_clock.tic();
        G->assemble_resistance();
        D->update<Statistics::AssembleVector>(t_clock.toc());

        if((initial_stiffness && counter != 0) || constant_matrix()) {
            t_clock.tic();
            if(SUANPAN_SUCCESS != G->process_load_resistance()) return SUANPAN_FAIL;
            if(SUANPAN_SUCCESS != G->process_constraint_resistance()) return SUANPAN_FAIL;
            D->update<Statistics::ProcessConstraint>(t_clock.toc());
        }
        else {
            t_clock.tic();
            G->assemble_matrix();
            D->update<Statistics::AssembleMatrix>(t_clock.toc());

            t_clock.tic();
            if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
            if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;
            D->update<Statistics::ProcessConstraint>(t_clock.toc());

            // indicate the global matrix has been assembled
            G->set_matrix_assembled_switch(true);
        }

        // call solver
        t_clock.tic();

        const auto residual = G->get_force_residual();

        auto flag = G->solve(samurai, residual);

        suanpan_assert([&] {
            if(!samurai.is_finite()) {
                suanpan_fatal("Non-finite number detected.\n");
                flag = SUANPAN_FAIL;
            }
        });

        // make sure system solver succeeds
        if(SUANPAN_SUCCESS != flag) {
            suanpan_error("Error code {} received.\n", flag);
            return flag;
        }

        inner_product = std::max(inner_product, std::fabs(dot(samurai, residual)));

        if(0u < counter)
            if(const auto amp = amplification(samurai, residual); amp > 0.) samurai *= amp;

        // deal with mpc
        if(const auto n_size = W->get_size(); 0 != W->get_multiplier_size()) {
            auto& aux_lambda = W->modify_auxiliary_lambda();
            auto& aux_border = W->get_auxiliary_stiffness();
            mat aux_right;
            if(SUANPAN_SUCCESS != G->solve(aux_right, aux_border)) return SUANPAN_FAIL;
            if(!solve(aux_lambda, aux_border.t() * aux_right.head_rows(n_size), aux_border.t() * samurai.head(n_size) - G->get_auxiliary_residual())) return SUANPAN_FAIL;
            samurai -= aux_right * aux_lambda;
        }

        D->update<Statistics::SolveSystem>(t_clock.toc());

        if(initial_stiffness) {
            if(!aitken) {
                aitken = true;
                pre_samurai = samurai;
            }
            else {
                aitken = false;
                const vec diff_samurai = pre_samurai - samurai;
                samurai *= dot(pre_samurai, diff_samurai) / dot(diff_samurai, diff_samurai);
            }
        }

        // avoid machine error accumulation
        G->erase_machine_error(samurai);

        // exit if converged
        // call corrector if it exists
        if(C->is_converged(counter)) return G->sync_status(true);
        // exit if maximum iteration is hit
        if(++counter > max_iteration) return SUANPAN_FAIL;

        // update internal variable
        G->update_internal(samurai);
        // update trial status for factory
        G->update_from_ninja();
        // for tracking
        G->update_load();
        // for tracking
        G->update_constraint();

        // fast handling for linear elastic case
        // sync status using newly computed increment across elements and nodes
        // this may just call predictor or call corrector
        if(D->get_attribute(ModalAttribute::LinearSystem)) return G->sync_status(false);
    }
}

void Newton::print() {
    suanpan_info("A solver based on Newton-Raphson method{}", initial_stiffness ? " using initial stiffness for each substep.\n" : ".\n");
}

double AICN::amplification(const vec& x, const vec& r) const {
    const auto hessian_norm = dot(x, r);
    if(hessian_norm <= 0.) return 0.;

    static constexpr auto ratio = .9;
    static constexpr auto factor = 2. / ratio * (1. / ratio - 1.);

    return 2. / (std::sqrt(1. + 2. * std::max(l_est, factor / std::sqrt(inner_product)) * std::sqrt(hessian_norm)) + 1.);
}

AICN::AICN(const unsigned T, const double L)
    : Newton(T)
    , l_est(fabs(L)) {}

void AICN::print() {
    suanpan_info("A solver based on cubic Newton method. arXiv:2211.00140\n");
}
