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

#include "GroupNodeRecorder.h"
#include <Domain/DomainBase.h>
#include <Domain/Group.h>

void GroupNodeRecorder::update_tag(const shared_ptr<DomainBase>& D) {
    vector<uword> tag;

    for(auto& I : groups) if(D->find_group(static_cast<unsigned>(I))) for(auto& J : D->get_group(static_cast<unsigned>(I))->get_pool()) tag.emplace_back(J);

    set_object_tag(unique(uvec(tag)));

    access::rw(get_data_pool()).resize(get_object_tag().n_elem);
}

GroupNodeRecorder::GroupNodeRecorder(const unsigned T, uvec&& B, const OutputType L, const unsigned I, const bool R, const bool H)
    : NodeRecorder(T, {}, L, I, R, H)
    , groups(std::forward<uvec>(B)) {}

void GroupNodeRecorder::initialize(const shared_ptr<DomainBase>& D) {
    update_tag(D);

    NodeRecorder::initialize(D);
}

void GroupNodeRecorder::print() { suanpan_info("A Node Recorder based on groups.\n"); }
