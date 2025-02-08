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

#include "Optimization.h"
#include <Domain/DomainBase.h>
#include <Solver/Integrator/Integrator.h>

int Optimization::initialize() {
    if(SUANPAN_SUCCESS != Static::initialize()) return SUANPAN_FAIL;

    if(0 == get_domain().lock()->get_criterion()) {
        suanpan_warning("At least one valid criterion shall be defined.\n");
        return SUANPAN_FAIL;
    }

    return SUANPAN_SUCCESS;
}

int Optimization::analyze() {
    const auto& D = get_domain().lock();
    auto& G = get_integrator();

    auto num_increment = 0u;

    while(true) {
        if(++num_increment > get_max_substep()) {
            suanpan_error("The maximum iteration {} reached.\n", get_max_substep());
            return SUANPAN_FAIL;
        }

        const auto code = Static::analyze();

        if(SUANPAN_FAIL == code) return SUANPAN_FAIL;
        if(SUANPAN_EXIT == code) return SUANPAN_SUCCESS;

        // const auto s_tag = D->get_current_step_tag();

        G->clear_status();

        // D->set_current_step_tag(s_tag);

        if(SUANPAN_SUCCESS != D->soft_restart()) return SUANPAN_FAIL;
    }
}
