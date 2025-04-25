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

#include "MaxDisplacement.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

unique_ptr<Criterion> MaxDisplacement::get_copy() { return make_unique<MaxDisplacement>(*this); }

int MaxDisplacement::process(const shared_ptr<DomainBase>& D) {
    const auto& t_vec = D->get_node(node)->get_reordered_dof();

    if(dof > t_vec.n_elem) return SUANPAN_SUCCESS;

    const auto& t_disp = D->get_factory()->get_current_displacement();

    return t_disp(t_vec.at(dof - 1)) > limit ? SUANPAN_EXIT : SUANPAN_SUCCESS;
}
