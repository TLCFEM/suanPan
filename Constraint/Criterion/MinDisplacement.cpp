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

#include "MinDisplacement.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

unique_ptr<Criterion> MinDisplacement::unique_copy() { return std::make_unique<MinDisplacement>(*this); }

int MinDisplacement::process(const shared_ptr<DomainBase>& D) {
    const auto t_vec = get_dof(D);

    if(t_vec.empty()) return SUANPAN_SUCCESS;

    const auto& t_disp = D->get_factory()->get_current_displacement();

    return t_disp(t_vec[0]) < limit ? SUANPAN_EXIT : SUANPAN_SUCCESS;
}
