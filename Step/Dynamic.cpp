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

#include "Dynamic.h"
#include <Converger/AbsIncreDisp.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Load/GroupNodalDisplacement.h>
#include <Solver/BFGS.h>
#include <Solver/Integrator/LeeNewmarkBase.h>
#include <Solver/MPDC.h>
#include <Solver/Newton.h>
#include <Solver/Ramm.h>

Dynamic::Dynamic(const unsigned T, const double P, const IntegratorType AT)
    : Step(T, P)
    , analysis_type(AT) {}

int Dynamic::initialize() {
    configure_storage_scheme();

    factory->set_analysis_type(AnalysisType::DYNAMICS);

    const auto& t_domain = database.lock();

    if(SUANPAN_SUCCESS != t_domain->restart()) return SUANPAN_FAIL;

    // converger
    if(nullptr == tester) tester = make_shared<AbsIncreDisp>();
    tester->set_domain(t_domain);

    // integrator
    if(nullptr == modifier) modifier = make_shared<Newmark>();
    else if(IntegratorType::Implicit == analysis_type) {
        if(IntegratorType::Implicit != modifier->type()) {
            suanpan_error("An implicit integrator is required.\n");
            return SUANPAN_FAIL;
        }
    }
    else if(IntegratorType::Implicit == modifier->type()) {
        suanpan_error("An explicit integrator is required.\n");
        return SUANPAN_FAIL;
    }
    modifier->set_domain(t_domain);

    // solver
    // avoid arc length solver
    if(nullptr != solver) if(dynamic_cast<Ramm*>(solver.get())) solver = nullptr;
    // automatically enable displacement controlled solver
    if(nullptr == solver) {
        auto flag = false;
        for(const auto& I : t_domain->get_load_pool())
            if(I->if_displacement_control() && I->get_start_step() == get_tag()) {
                flag = true;
                break;
            }
        flag ? solver = make_shared<MPDC>() : solver = make_shared<Newton>();
    }

    if(dynamic_cast<BFGS*>(solver.get()) && dynamic_cast<LeeNewmarkBase*>(modifier.get())) {
        suanpan_error("BFGS solver is not supported by Lee's damping model.\n");
        return SUANPAN_FAIL;
    }

    solver->set_converger(tester);
    solver->set_integrator(modifier);

    if(SUANPAN_SUCCESS != tester->initialize()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != modifier->initialize()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != solver->initialize()) return SUANPAN_FAIL;

    return SUANPAN_SUCCESS;
}

int Dynamic::analyze() {
    auto& S = get_solver();
    auto& G = get_integrator();
    auto& W = get_factory();

    auto remain_time = get_time_period();
    auto step_time = get_ini_step_size();

    unsigned num_increment = 0, num_converged_step = 0;

    while(true) {
        // check if the target time point is hit
        if(remain_time <= 1E-7) return SUANPAN_SUCCESS;
        // check if the maximum substep number is hit
        if(++num_increment > get_max_substep()) {
            suanpan_warning("The maximum sub-step number {} reached.\n", get_max_substep());
            return SUANPAN_FAIL;
        }
        // update incremental and trial time
        G->update_incre_time(step_time);
        if(const auto code = S->analyze(); SUANPAN_SUCCESS == code) {
            // success step
            // eat current increment
            set_time_left(remain_time -= W->get_incre_time());
            // commit converged iteration
            G->stage_and_commit_status();
            // record response
            G->record();
            if(G->allow_to_change_time_step()) {
                if(!is_fixed_step_size() && ++num_converged_step > 5) {
                    step_time = std::min(get_max_step_size(), step_time * time_step_amplification);
                    num_converged_step = 0;
                }
                // check if time overflows
                if(step_time > remain_time) step_time = remain_time;
            }
        }
        else if(SUANPAN_FAIL == code) {
            // failed step
            // reset to the start of current substep
            G->reset_status();
            // check if minimum step size is hit
            if(step_time <= get_min_step_size()) {
                suanpan_error("The minimum step size {:.3E} reached.\n", get_min_step_size());
                return SUANPAN_FAIL;
            }
            // check if fixed step size
            if(is_fixed_step_size() || !G->allow_to_change_time_step()) {
                suanpan_error("Cannot converge with the given step size {:.3E}.\n", step_time);
                return SUANPAN_FAIL;
            }
            // step size is allowed to decrease
            step_time *= .5;
            set_max_substep(num_increment + static_cast<unsigned>(remain_time / step_time) + 1);
            num_converged_step = 0;
        }
        else return SUANPAN_FAIL; // positive codes are from lapack subroutines
    }
}
