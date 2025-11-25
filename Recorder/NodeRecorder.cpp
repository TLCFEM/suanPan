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

#include "NodeRecorder.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

void NodeRecorder::record_impl(const shared_ptr<DomainBase>& D) {
    const auto populate = [&](const vec& in) {
        for(const auto I : object_tag) insert({in(D->get<Node>(I)->get_reordered_dof())}, I);
    };

    if(OutputType::GDF == variable_type) {
        auto& damping_force = D->get_factory()->get_current_damping_force();
        if(damping_force.empty()) return;
        populate(damping_force);
    }
    else if(OutputType::GIF == variable_type) {
        auto& inertial_force = D->get_factory()->get_current_inertial_force();
        if(inertial_force.empty()) return;
        populate(inertial_force);
    }
    else if(OutputType::MM == variable_type) {
        auto& momentum = D->get_factory()->get_momentum();
        if(momentum.empty()) return;
        populate(momentum);
    }
    else
        for(const auto I : object_tag) insert(D->get<Node>(I)->record(variable_type), I);

    insert(D->get_factory()->get_current_time());
}

void NodeRecorder::initialize(const shared_ptr<DomainBase>& D) {
    update_tag(D);

    std::vector<uword> pool;
    pool.reserve(object_tag.n_elem);
    for(const auto I : object_tag)
        if(!D->find<Node>(I) || !D->get<Node>(I)->is_active())
            suanpan_warning("Node {} is not available/active, removed from recorder {}.\n", I, get_tag());
        else pool.emplace_back(I);

    object_tag = pool;
}

void NodeRecorder::print() { suanpan_info("A node recorder.\n"); }

const uvec& GroupNodeRecorder::update_tag(const shared_ptr<DomainBase>& D) { return object_tag = D->flatten_group(reference_tag); }

void GroupNodeRecorder::print() { suanpan_info("A node recorder based on groups.\n"); }
