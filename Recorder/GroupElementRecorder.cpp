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

#include "GroupElementRecorder.h"
#include <Domain/DomainBase.h>
#include <Domain/Group/Group.h>

void GroupElementRecorder::update_tag(const shared_ptr<DomainBase>& D) {
    std::vector<uword> tag;

    for(const auto I : groups) if(D->find<Group>(I)) for(const auto J : D->get<Group>(I)->get_pool()) tag.emplace_back(J);

    set_object_tag(unique(uvec(tag)));

    access::rw(get_data_pool()).resize(get_object_tag().n_elem);
}

GroupElementRecorder::GroupElementRecorder(const unsigned T, uvec&& B, const OutputType L, const unsigned I, const bool R, const bool H)
    : ElementRecorder(T, {}, L, I, R, H)
    , groups(std::move(B)) {}

void GroupElementRecorder::initialize(const shared_ptr<DomainBase>& D) {
    update_tag(D);

    ElementRecorder::initialize(D);
}

void GroupElementRecorder::print() {
    suanpan_info("An element recorder based on groups.\n");
}
