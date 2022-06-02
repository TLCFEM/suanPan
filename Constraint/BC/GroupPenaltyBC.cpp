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

#include "GroupPenaltyBC.h"
#include <Domain/DomainBase.h>
#include <Domain/Group/Group.h>

void GroupPenaltyBC::update_node_tag(const shared_ptr<DomainBase>& D) {
    vector<uword> tag;

    for(const auto I : groups) if(D->find<Group>(I)) for(const auto J : D->get<Group>(I)->get_pool()) tag.emplace_back(J);

    node_encoding = unique(uvec(tag));
}

GroupPenaltyBC::GroupPenaltyBC(const unsigned T, const unsigned S, uvec&& N, uvec&& D)
    : MultiplierBC(T, S, uvec{}, std::forward<uvec>(D))
    , groups(std::forward<uvec>(N)) {}

GroupPenaltyBC::GroupPenaltyBC(const unsigned T, const unsigned S, uvec&& N, const char TP)
    : MultiplierBC(T, S, uvec{}, TP)
    , groups(std::forward<uvec>(N)) {}

int GroupPenaltyBC::initialize(const shared_ptr<DomainBase>& D) {
    update_node_tag(D);

    return MultiplierBC::initialize(D);
}

int GroupPenaltyBC::process(const shared_ptr<DomainBase>& D) {
    return PenaltyBC::process(D); // NOLINT(bugprone-parent-virtual-call)
}

int GroupPenaltyBC::process_resistance(const shared_ptr<DomainBase>& D) {
    return PenaltyBC::process_resistance(D); // NOLINT(bugprone-parent-virtual-call)
}
