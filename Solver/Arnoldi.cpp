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

#include "Arnoldi.h"
#include <Domain/DomainBase.h>
#include <Domain/FactoryHelper.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Toolbox/arpack_wrapper.h>

Arnoldi::Arnoldi(const unsigned T, const unsigned N, const char TP)
    : Solver(T)
    , eigen_num(N)
    , eigen_type(TP) {}

int Arnoldi::initialize() {
    if(get_integrator() == nullptr) {
        suanpan_error("initialize() needs a valid integrator.\n");
        return SUANPAN_FAIL;
    }

    return SUANPAN_SUCCESS;
}

int Arnoldi::analyze() {
    auto& G = get_integrator();
    const auto& D = G->get_domain();
    auto& W = D->get_factory();

    if(SUANPAN_SUCCESS != G->process_modifier()) return SUANPAN_FAIL;

    D->assemble_trial_mass();
    D->assemble_trial_stiffness();

    // if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;

    return eig_solve(get_eigenvalue(W), get_eigenvector(W), W->get_stiffness(), W->get_mass(), eigen_num, 'L' == eigen_type ? "LM" : "SM");
}

void Arnoldi::print() { suanpan_info("A solver using Arnoldi method.\n"); }
