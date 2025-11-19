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

void NodeRecorder::initialize(const shared_ptr<DomainBase>& D) {
    std::vector<uword> pool;
    pool.reserve(get_object_tag().n_elem);
    for(const auto I : get_object_tag())
        if(!D->find<Node>(I) || !D->get<Node>(I)->is_active())
            suanpan_warning("Node {} is not available/active, removed from recorder {}.\n", I, get_tag());
        else pool.emplace_back(I);

    set_object_tag(pool);

    access::rw(get_data_pool()).resize(get_object_tag().n_elem);
}

void NodeRecorder::record(const shared_ptr<DomainBase>& D) {
    if(!if_perform_record()) return;

    auto& obj_tag = get_object_tag();

    if(OutputType::GDF == get_variable_type()) {
        auto& damping_force = D->get_factory()->get_current_damping_force();
        if(damping_force.empty()) return;

        for(unsigned I = 0; I < obj_tag.n_elem; ++I)
            if(const auto& t_node = D->get<Node>(obj_tag(I)); t_node->is_active()) insert({damping_force(t_node->get_reordered_dof())}, I);
    }
    else if(OutputType::GIF == get_variable_type()) {
        auto& inertial_force = D->get_factory()->get_current_inertial_force();
        if(inertial_force.empty()) return;

        for(unsigned I = 0; I < obj_tag.n_elem; ++I)
            if(const auto& t_node = D->get<Node>(obj_tag(I)); t_node->is_active()) insert({inertial_force(t_node->get_reordered_dof())}, I);
    }
    else if(OutputType::MM == get_variable_type()) {
        auto& momentum = D->get_factory()->get_momentum();
        if(momentum.empty()) return;

        for(unsigned I = 0; I < obj_tag.n_elem; ++I)
            if(const auto& t_node = D->get<Node>(obj_tag(I)); t_node->is_active()) insert({momentum(t_node->get_reordered_dof())}, I);
    }
    else
        for(unsigned I = 0; I < obj_tag.n_elem; ++I)
            if(const auto& t_node = D->get<Node>(obj_tag(I)); t_node->is_active()) insert(t_node->record(get_variable_type()), I);

    if(if_record_time()) insert(D->get_factory()->get_current_time());
}

void NodeRecorder::print() {
    suanpan_info("A node recorder.\n");
}
