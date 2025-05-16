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

SumRecorder::SumRecorder(const unsigned T, uvec&& B, const OutputType L, const unsigned I, const bool R, const bool H)
    : Recorder(T, std::move(B), L, I, R, H) { access::rw(get_data_pool()).resize(1); }

void SumRecorder::initialize(const shared_ptr<DomainBase>& D) {
    for(const auto I : get_object_tag())
        if(!D->find<Node>(I)) {
            D->disable_recorder(get_tag());
            return;
        }
}

void SumRecorder::record(const shared_ptr<DomainBase>& D) {
    if(!if_perform_record()) return;

    auto data = 0.;
    for(const auto I : get_object_tag()) {
        const auto& t_node = D->get<Node>(I);
        if(!t_node->is_active()) continue;
        const auto t_data = t_node->record(get_variable_type());
        if(t_data.empty()) continue;
        const auto& n_data = t_data[0];
        if(n_data.empty()) continue;
        data += as_scalar(n_data);
    }
    insert({{data}}, 0);

    if(if_record_time()) insert(D->get_factory()->get_current_time());
}

void SumRecorder::print() {
    suanpan_info("A summation recorder computes the summation of a collection of nodal scalar variables.\n");
}
