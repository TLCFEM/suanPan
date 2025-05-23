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

#include "GroupGroup.h"

#include <Domain/DomainBase.h>
#include <Toolbox/utility.h>

GroupGroup::GroupGroup(const unsigned T, uvec&& GT)
    : Group(T)
    , group_tag(std::move(GT)) {}

void GroupGroup::initialize(const shared_ptr<DomainBase>& D) {
    uword size = 0;
    std::vector<const uvec*> ocean;
    ocean.reserve(group_tag.n_elem);
    for(const auto I : group_tag)
        if(D->find<Group>(I)) {
            auto& t_group = D->get<Group>(I);
            t_group->initialize(D);
            ocean.emplace_back(&t_group->get_pool());
            size += ocean.back()->size();
        }

    std::vector<uword> pond;
    pond.reserve(size);

    for(const auto& I : ocean)
        for(auto J : *I) pond.emplace_back(J);

    pool = suanpan::unique(pond);
}
