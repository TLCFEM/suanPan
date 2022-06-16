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

#include "Group.h"

Group::Group(const unsigned T)
    : Tag(T) { suanpan_debug("Group %u ctor() called.\n", get_tag()); }

Group::Group(const unsigned T, uvec&& R)
    : Tag(T)
    , pool(std::forward<uvec>(R)) { suanpan_debug("Group %u ctor() called.\n", get_tag()); }

Group::~Group() { suanpan_debug("Group %u dtor() called.\n", get_tag()); }

void Group::initialize(const shared_ptr<DomainBase>&) {}

const uvec& Group::get_pool() const { return pool; }

void Group::print() {
    suanpan_info("A Group object with tag %u contains the following tags:\n", get_tag());
    pool.t().print();
}
