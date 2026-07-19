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

#include "MinResistance.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

unique_ptr<Criterion> MinResistance::unique_copy() { return std::make_unique<MinResistance>(*this); }

int MinResistance::process(const shared_ptr<DomainBase>& D) {
    if(const auto t_vec = get_dof(D); !t_vec.empty()) return D->get_factory()->get_current_resistance()(t_vec[0]) < limit ? SUANPAN_EXIT : SUANPAN_SUCCESS;

    return SUANPAN_SUCCESS;
}
