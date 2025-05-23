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

#include "Buckle.h"

#include <Converger/AbsIncreDisp.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Newton.h>
#include <Toolbox/arpack.h>

Buckle::Buckle(const unsigned T)
    : Static(T, 1.) {}

int Buckle::initialize() {
#ifdef SUANPAN_DISTRIBUTED
    suanpan_warning("Buckling analysis currently does not support distributed computation thus it will be conducted on each node.\n");
#endif

    const auto t_domain = database.lock();

    // converger
    if(nullptr == tester) tester = std::make_shared<AbsIncreDisp>();
    tester->set_domain(t_domain);

    // integrator
    modifier = std::make_shared<Integrator>();
    modifier->set_domain(t_domain);

    // solver
    if(nullptr == solver) solver = std::make_shared<Newton>();
    solver->set_converger(tester);
    solver->set_integrator(modifier);

    if(SUANPAN_SUCCESS != tester->initialize()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != modifier->initialize()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != solver->initialize()) return SUANPAN_FAIL;

    configure_storage_scheme();

    factory->set_analysis_type(AnalysisType::BUCKLE);

    return t_domain->restart();
}

int Buckle::analyze() {
    if(Static::analyze() == SUANPAN_FAIL) return SUANPAN_FAIL;

    const auto D = get_domain().lock();
    auto& G = get_integrator();
    auto& W = get_factory();

    if(!W->is_nlgeom()) {
        suanpan_error("Buckling analysis requires an active nlgeom model.\n");
        return SUANPAN_FAIL;
    }

    // assemble stiffness and geometry as they may be modified in solver
    D->assemble_trial_stiffness();
    D->assemble_trial_geometry();

    if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;

    const shared_ptr t_geometry = W->get_geometry()->make_copy();
    t_geometry *= -1.;

    if(eig_solve(W->modify_eigenvalue(), W->modify_eigenvector(), W->get_stiffness(), t_geometry, 1, "SM") != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    suanpan_info("\nBuckling load multiplier: {:.8E}.\n", W->get_eigenvalue().at(0));

    // record response
    G->record();

    return SUANPAN_SUCCESS;
}
