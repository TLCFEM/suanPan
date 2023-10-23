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

#include "ArcLength.h"
#include <Converger/AbsIncreDisp.h>
#include <Domain/Domain.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Ramm.h>

ArcLength::ArcLength(const unsigned T, const unsigned NT, const unsigned DT, const double MA)
    : Step(T, 0.)
    , node(NT)
    , dof(DT)
    , magnitude(MA) {}

int ArcLength::initialize() {
    const auto t_domain = database.lock();

    // converger
    if(nullptr == tester) tester = make_shared<AbsIncreDisp>();
    tester->set_domain(t_domain);

    // integrator
    modifier = make_shared<Integrator>();
    modifier->set_domain(t_domain);

    // solver
    if(nullptr != solver && nullptr == std::dynamic_pointer_cast<Ramm>(solver)) solver = nullptr;
    if(nullptr == solver) solver = make_shared<Ramm>();
    solver->set_converger(tester);
    solver->set_integrator(modifier);

    if(SUANPAN_SUCCESS != tester->initialize()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != modifier->initialize()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != solver->initialize()) return SUANPAN_FAIL;

    configure_storage_scheme();

    factory->set_solver_type(sparse_mat ? SolverType::MUMPS : SolverType::LAPACK);

    factory->set_analysis_type(AnalysisType::STATICS);

    if(SUANPAN_SUCCESS != t_domain->restart()) return SUANPAN_FAIL;

    factory->set_reference_size(1);
    factory->initialize_load_factor();

    factory->modify_reference_load()(t_domain->get_node(node)->get_reordered_dof().at(dof - 1)) = magnitude;

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
