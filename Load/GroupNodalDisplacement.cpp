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

#include "GroupNodalDisplacement.h"
#include <Domain/DomainBase.h>
#include <Domain/Group.h>

void GroupNodalDisplacement::update_node_tag(const shared_ptr<DomainBase>& D) {
    vector<uword> tag;

    for(auto& I : groups) if(D->find<Group>(I)) for(auto& J : D->get<Group>(I)->get_pool()) tag.emplace_back(J);

    node_encoding = unique(uvec(tag));
}

GroupNodalDisplacement::GroupNodalDisplacement(const unsigned T, const unsigned ST, const double L, uvec&& N, const unsigned D, const unsigned AT)
    : NodalDisplacement(T, ST, L, uvec{}, uvec{D}, AT)
    , groups(std::forward<uvec>(N)) {}

GroupNodalDisplacement::GroupNodalDisplacement(const unsigned T, const unsigned ST, const double L, uvec&& N, uvec&& D, const unsigned AT)
    : NodalDisplacement(T, ST, L, uvec{}, std::forward<uvec>(D), AT)
    , groups(std::forward<uvec>(N)) {}

int GroupNodalDisplacement::initialize(const shared_ptr<DomainBase>& D) {
    update_node_tag(D);

    return NodalDisplacement::initialize(D);
}
