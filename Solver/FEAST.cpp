/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#include "FEAST.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Toolbox/eigen_solver.h>

FEAST::FEAST(const unsigned T, const unsigned N, const double C, const double R)
    : Solver(T)
    , eigen_num(N)
    , centre(C)
    , radius(R) {}

int FEAST::initialize() {
    auto& G = get_integrator();

    if(nullptr == G) {
        suanpan_error("A valid integrator is required.\n");
        return SUANPAN_FAIL;
    }

    const auto scheme = G->get_domain()->get_factory()->get_storage_scheme();

#ifdef SUANPAN_DISTRIBUTED
    if(StorageScheme::SPARSE != scheme && StorageScheme::SPARSESYMM != scheme) {
        suanpan_error("The parallel FEAST solver supports sparse storage only.\n");
        return SUANPAN_FAIL;
    }
#else
    if(StorageScheme::SYMMPACK == scheme) {
        suanpan_error("The symmetric pack storage is not supported.\n");
        return SUANPAN_FAIL;
    }
#endif

    return SUANPAN_SUCCESS;
}

int FEAST::analyze() {
    auto& G = get_integrator();
    const auto D = G->get_domain();
    auto& W = D->get_factory();

    wall_clock t_clock;

    t_clock.tic();
    if(SUANPAN_SUCCESS != G->process_modifier()) return SUANPAN_FAIL;
    D->update<Statistics::UpdateStatus>(t_clock);

    t_clock.tic();
    D->assemble_trial_mass();
    D->assemble_trial_stiffness();
    D->update<Statistics::AssembleMatrix>(t_clock);

    t_clock.tic();
    // if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;
    D->update<Statistics::ProcessConstraint>(t_clock);

    return feast_solve(W, eigen_num, centre, radius);
}

void FEAST::print() { suanpan_info("An eigen solver using FEAST solver.\n"); }
