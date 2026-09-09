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

#include "Criterion.h"

#include <Domain/DomainBase.h>

void Criterion::set_start_step(const unsigned ST) { start_step = ST; }

void Criterion::set_end_step(const unsigned ST) { end_step = ST; }

bool Criterion::if_apply(const shared_ptr<DomainBase>& D) const {
    const auto t_step = D->get_current_step_tag();
    return start_step != 0u && t_step >= start_step && t_step < end_step && is_active();
}

int Criterion::initialize(const shared_ptr<DomainBase>&) { return SUANPAN_SUCCESS; }
