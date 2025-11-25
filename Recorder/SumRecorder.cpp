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

#include "SumRecorder.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

void SumRecorder::record_impl(const shared_ptr<DomainBase>& D) {
    running_stat_vec<vec> data;
    for(const auto I : object_tag)
        for(const auto& t_data : D->get<Node>(I)->record(variable_type)) data(t_data);

    insert({data.mean() * data.count()}, 0);

    insert(D->get_factory()->get_current_time());
}

void SumRecorder::initialize(const shared_ptr<DomainBase>& D) {
    for(const auto I : update_tag(D))
        if(!D->find<Node>(I) || !D->get<Node>(I)->is_active()) {
            suanpan_warning("Node {} is not available/active, recorder {} is disabled.\n", I, get_tag());
            D->disable_recorder(get_tag());
            return;
        }
}

void SumRecorder::print() { suanpan_info("A summation recorder computes the summation of a collection of nodal scalar variables.\n"); }

const uvec& GroupSumRecorder::update_tag(const shared_ptr<DomainBase>& D) { return object_tag = D->flatten_group(reference_tag); }

void GroupSumRecorder::print() { suanpan_info("A summation recorder.\n"); }
