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

#include "Load.h"

#include <Domain/DomainBase.h>
#include <Domain/Group/Group.h>

double Load::multiplier = 1E8;

Load::Load(const unsigned T, const unsigned ST, const unsigned AT, uvec&& NT, uvec&& DT, const double PT)
    : ConditionalModifier(T, AT, std::move(NT), std::move(DT))
    , pattern(PT) {}

void Load::enable_displacement_control() const { access::rw(mpdc_flag) = true; }

bool Load::if_displacement_control() const { return mpdc_flag; }

const vec& Load::get_trial_load() const { return trial_load; }

const vec& Load::get_trial_settlement() const { return trial_settlement; }

const sp_vec& Load::get_reference_load() const { return reference_load; }

void set_load_multiplier(const double M) { Load::multiplier = M; }

GroupLoad::GroupLoad(uvec&& N)
    : groups(std::move(N)) {}

uvec GroupLoad::update_object_tag(const shared_ptr<DomainBase>& D) const {
    suanpan::unordered_set<uword> tag;

    for(const auto I : groups) {
        const auto& t_group = D->get<Group>(I);
        if(nullptr == t_group) continue;
        const auto& t_pool = t_group->get_pool();
        tag.insert(t_pool.cbegin(), t_pool.cend());
    }

    return to_uvec(tag);
}
