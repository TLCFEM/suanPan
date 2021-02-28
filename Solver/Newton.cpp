////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2021 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

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
	auto& W = G->get_domain().lock()->get_factory();

	suanpan_info("current analysis time: %.5f.\n", W->get_trial_time());

	const auto& max_iteration = C->get_max_iteration();

	auto flag = 0;

	// iteration counter
	unsigned counter = 0;

	// ninja alias
	auto& ninja = get_ninja(W);

	vec pre_ninja;

	auto aitken = false;

	while(true) {
		// process modifiers
		if(G->process_modifier() != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		// assemble resistance
		G->assemble_resistance();

		if(initial_stiffness && counter != 0) flag = G->solve_trs(ninja, G->get_force_residual());
		else {
			// first iteration
			// assemble stiffness
			G->assemble_matrix();
			// process loads
			if(G->process_load() != SUANPAN_SUCCESS) return SUANPAN_FAIL;
			// process constraints
			if(G->process_constraint() != SUANPAN_SUCCESS) return SUANPAN_FAIL;
			// call solver
			flag = G->solve(ninja, G->get_force_residual());
		}

		suanpan_debug([&]() {
			if(!ninja.is_finite()) {
				suanpan_fatal("infinite number detected.\n");
				flag = SUANPAN_FAIL;
			}
		});

		// make sure lapack solver succeeds
		if(SUANPAN_SUCCESS != flag) return flag;

		// deal with mpc
		if(const auto n_size = W->get_size(); 0 != W->get_mpc()) {
			auto& border = W->get_auxiliary_stiffness();
			mat right;
			if(G->solve_trs(right, border) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
			vec aux_factor;
			if(!solve(aux_factor, border.t() * right.head_rows(n_size), border.t() * ninja.head_rows(n_size) - W->get_auxiliary_load())) return SUANPAN_FAIL;
			// W->update_trial_auxiliary_resistance(W->get_trial_auxiliary_resistance() + border * aux_factor);
			ninja -= right * aux_factor;
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
		// update internal variable
		G->update_internal(ninja);
		// update trial status for factory
		W->update_trial_displacement(W->get_trial_displacement() + ninja);
		// update for nodes and elements
		if(G->update_trial_status() != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		// exit if converged
		if(C->is_converged()) return SUANPAN_SUCCESS;
		// exit if maximum iteration is hit
		if(++counter > max_iteration) return SUANPAN_FAIL;
	}
}

void Newton::print() { suanpan_info("A solver based on Newton--Raphson iteration method %s", initial_stiffness ? "using initial stiffness for each substep.\n" : ".\n"); }
