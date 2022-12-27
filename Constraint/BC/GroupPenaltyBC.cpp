/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

GroupPenaltyBC::GroupPenaltyBC(const unsigned T, const unsigned S, uvec&& N, uvec&& D)
    : MultiplierBC(T, S, uvec{}, std::forward<uvec>(D))
    , groups(std::forward<uvec>(N)) {}

GroupPenaltyBC::GroupPenaltyBC(const unsigned T, const unsigned S, uvec&& N, const char TP)
    : MultiplierBC(T, S, uvec{}, TP)
    , groups(std::forward<uvec>(N)) {}

int GroupPenaltyBC::initialize(const shared_ptr<DomainBase>& D) {
    suanpan::unordered_set<uword> tag;

    for(const auto I : groups) {
        const auto& t_group = D->get<Group>(I);
        if(nullptr == t_group) continue;
        const auto& t_pool = t_group->get_pool();
        tag.insert(t_pool.cbegin(), t_pool.cend());
    }

    node_encoding = to_uvec(tag);

    return MultiplierBC::initialize(D);
}

int GroupPenaltyBC::process(const shared_ptr<DomainBase>& D) {
    return PenaltyBC::process(D); // NOLINT(bugprone-parent-virtual-call)
}

int GroupPenaltyBC::process_resistance(const shared_ptr<DomainBase>& D) {
    return PenaltyBC::process_resistance(D); // NOLINT(bugprone-parent-virtual-call)
}
