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

#include "Criterion.h"

#include <Domain/DomainBase.h>

Criterion::Criterion(const unsigned T, const unsigned ST)
    : CopiableTag(T)
    , step_tag(ST) {}

void Criterion::set_step_tag(const unsigned T) { step_tag = T; }

unsigned Criterion::get_step_tag() const { return step_tag; }

int Criterion::initialize(const shared_ptr<DomainBase>&) { return SUANPAN_SUCCESS; }
