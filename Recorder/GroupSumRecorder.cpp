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

#include "GroupSumRecorder.h"
#include <Domain/DomainBase.h>
#include <Domain/Group/Group.h>

void GroupSumRecorder::update_tag(const shared_ptr<DomainBase>& D) {
    std::vector<uword> tag;

    for(const auto I : groups) if(D->find<Group>(I)) for(const auto J : D->get<Group>(I)->get_pool()) tag.emplace_back(J);

    set_object_tag(unique(uvec(tag)));
}

GroupSumRecorder::GroupSumRecorder(const unsigned T, uvec&& B, const OutputType L, const unsigned I, const bool R, const bool H)
    : SumRecorder(T, {}, L, I, R, H)
    , groups(std::move(B)) { access::rw(get_data_pool()).resize(1); }

void GroupSumRecorder::initialize(const shared_ptr<DomainBase>& D) {
    update_tag(D);

    SumRecorder::initialize(D);
}

void GroupSumRecorder::print() {
    suanpan_info("A summation recorder.\n");
}
