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

#include "Solver.h"
#include <Converger/Converger.h>
#include <Solver/Integrator/Integrator.h>
#include <Domain/DomainBase.h>
#include <Step/Step.h>

Solver::Solver(const unsigned T)
    : Tag(T) { suanpan_debug("Solver %u ctor() called.\n", get_tag()); }

Solver::~Solver() { suanpan_debug("Solver %u dtor() called.\n", get_tag()); }

int Solver::initialize() {
    if(nullptr == converger) {
        suanpan_error("initialize() needs a valid converger.\n");
        return SUANPAN_FAIL;
    }

    if(nullptr == modifier) {
        suanpan_error("initialize() needs a valid integrator.\n");
        return SUANPAN_FAIL;
    }

    return SUANPAN_SUCCESS;
}

void Solver::set_converger(const shared_ptr<Converger>& C) { converger = C; }

const shared_ptr<Converger>& Solver::get_converger() const { return converger; }

void Solver::set_integrator(const shared_ptr<Integrator>& G) { modifier = G; }

const shared_ptr<Integrator>& Solver::get_integrator() const { return modifier; }

bool Solver::constant_matrix() const {
    auto& G = get_integrator();
    const auto& D = G->get_domain();
    auto& S = D->get_current_step();

    // need to satisfy a number of conditions:
    // 1. fixed step size
    // 2. if not fixed step size, the effective stiffness needs to be independent from time
    // 3. the system needs to be linear
    // 4. the effective stiffness has been assembled
    return (S->is_fixed_step_size() || G->time_independent_matrix()) && D->get_attribute(ModalAttribute::LinearSystem) && G->matrix_is_assembled();
}
