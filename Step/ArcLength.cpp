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

#include "ArcLength.h"

#include <Converger/AbsIncreDisp.h>
#include <Domain/Domain.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Ramm.h>

ArcLength::ArcLength(const unsigned T)
    : Step(T, 0.) {}

int ArcLength::initialize() {
#ifdef SUANPAN_DISTRIBUTED
    suanpan_error("Arc-length analysis currently does not support distributed computation.\n");
    return SUANPAN_FAIL;
#endif

    const auto t_domain = database.lock();

    // converger
    if(nullptr == tester) tester = std::make_shared<AbsIncreDisp>();
    tester->set_domain(t_domain);

    // integrator
    modifier = std::make_shared<Integrator>();
    modifier->set_domain(t_domain);

    // solver
    if(nullptr == solver || nullptr == std::dynamic_pointer_cast<Ramm>(solver)) solver = std::make_shared<Ramm>();
    solver->set_converger(tester);
    solver->set_integrator(modifier);

    if(SUANPAN_SUCCESS != tester->initialize()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != modifier->initialize()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != solver->initialize()) return SUANPAN_FAIL;

    configure_storage_scheme();

    factory->set_solver_type(sparse_mat ? SolverType::SUPERLU : SolverType::LAPACK);

    factory->set_analysis_type(AnalysisType::STATICS);

    if(SUANPAN_SUCCESS != t_domain->restart()) return SUANPAN_FAIL;

    factory->set_reference_size(1);
    factory->initialize_load_factor();

    return SUANPAN_SUCCESS;
}

int ArcLength::analyze() {
    auto& S = get_solver();
    auto& G = get_integrator();

    auto num_iteration = 0u;

    auto arc_length = get_ini_step_size();
    const auto min_arc_length = get_min_step_size();
    const auto max_arc_length = get_max_step_size();

    while(true) {
        if(num_iteration++ > get_max_substep()) {
            suanpan_warning("The maximum sub-step number {} reached.\n", get_max_substep());
            return SUANPAN_FAIL;
        }
        S->set_step_size(arc_length);
        if(auto code = S->analyze(); code == SUANPAN_SUCCESS) {
            G->stage_and_commit_status();
            G->record();
            // adjust arc length, always increase
            if(!is_fixed_step_size()) {
                arc_length *= S->get_step_amplifier();
                if(max_arc_length > 0. && arc_length > max_arc_length) arc_length = max_arc_length;
            }
            // if exit is returned, the analysis shall be terminated
            code = G->process_criterion();
            if(SUANPAN_SUCCESS != code) return code;
        }
        else if(code == SUANPAN_FAIL) {
            G->reset_status();
            // no way to converge
            if(is_fixed_step_size() || suanpan::approx_equal(arc_length, min_arc_length)) return SUANPAN_FAIL;
            // decrease arc length
            arc_length *= .5;
            if(min_arc_length > 0. && arc_length < min_arc_length) arc_length = min_arc_length;
        }
        else return SUANPAN_FAIL;
    }
}
