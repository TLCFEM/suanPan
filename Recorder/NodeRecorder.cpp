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

#include <Domain/DOF.h>
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

    auto insert_damping_force = [&](const uword J) {
        for(unsigned I = 0; I < obj_tag.n_elem; ++I) {
            const auto& t_node = D->get<Node>(obj_tag(I));
            if(!t_node->is_active()) continue;
            const auto& t_dof = t_node->get_reordered_dof();
            const auto& t_force = D->get_factory()->get_current_damping_force();
            insert({{(t_dof.n_elem > J && t_force.n_elem > t_dof(J) ? t_force(t_dof(J)) : 0.)}}, I);
        }
    };
    auto insert_inertial_force = [&](const uword J) {
        for(unsigned I = 0; I < obj_tag.n_elem; ++I) {
            const auto& t_node = D->get<Node>(obj_tag(I));
            if(!t_node->is_active()) continue;
            const auto& t_dof = t_node->get_reordered_dof();
            const auto& t_force = D->get_factory()->get_current_inertial_force();
            insert({{(t_dof.n_elem > J && t_force.n_elem > t_dof(J) ? t_force(t_dof(J)) : 0.)}}, I);
        }
    };

    auto insert_momentum = [&](const DOF DI) {
        for(unsigned I = 0; I < obj_tag.n_elem; ++I) {
            const auto& t_node = D->get<Node>(obj_tag(I));
            if(!t_node->is_active()) continue;
            const auto& t_dof = t_node->get_reordered_dof();
            const auto& t_dof_identifier = t_node->get_dof_identifier();
            const auto& t_momentum = D->get_factory()->get_momentum();
            const auto [flag, position] = if_contain(t_dof_identifier, DI);
            insert({{(flag && t_momentum.n_elem > t_dof(position) ? t_momentum(t_dof(position)) : 0.)}}, I);
        }
    };

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
    else if(OutputType::GDF1 == get_variable_type()) insert_damping_force(0);
    else if(OutputType::GDF2 == get_variable_type()) insert_damping_force(1);
    else if(OutputType::GDF3 == get_variable_type()) insert_damping_force(2);
    else if(OutputType::GDF4 == get_variable_type() || OutputType::GDM1 == get_variable_type()) insert_damping_force(3);
    else if(OutputType::GDF5 == get_variable_type() || OutputType::GDM2 == get_variable_type()) insert_damping_force(4);
    else if(OutputType::GDF6 == get_variable_type() || OutputType::GDM3 == get_variable_type()) insert_damping_force(5);
    else if(OutputType::GIF1 == get_variable_type()) insert_inertial_force(0);
    else if(OutputType::GIF2 == get_variable_type()) insert_inertial_force(1);
    else if(OutputType::GIF3 == get_variable_type()) insert_inertial_force(2);
    else if(OutputType::GIF4 == get_variable_type() || OutputType::GIM1 == get_variable_type()) insert_inertial_force(3);
    else if(OutputType::GIF5 == get_variable_type() || OutputType::GIM2 == get_variable_type()) insert_inertial_force(4);
    else if(OutputType::GIF6 == get_variable_type() || OutputType::GIM3 == get_variable_type()) insert_inertial_force(5);
    else if(OutputType::MM1 == get_variable_type()) insert_momentum(DOF::U1);
    else if(OutputType::MM2 == get_variable_type()) insert_momentum(DOF::U2);
    else if(OutputType::MM3 == get_variable_type()) insert_momentum(DOF::U3);
    else if(OutputType::MM4 == get_variable_type() || OutputType::MMR1 == get_variable_type()) insert_momentum(DOF::UR1);
    else if(OutputType::MM5 == get_variable_type() || OutputType::MMR2 == get_variable_type()) insert_momentum(DOF::UR2);
    else if(OutputType::MM6 == get_variable_type() || OutputType::MMR3 == get_variable_type()) insert_momentum(DOF::UR3);
    else
        for(unsigned I = 0; I < obj_tag.n_elem; ++I)
            if(const auto& t_node = D->get<Node>(obj_tag(I)); t_node->is_active()) insert(t_node->record(get_variable_type()), I);

    if(if_record_time()) insert(D->get_factory()->get_current_time());
}

void NodeRecorder::print() {
    suanpan_info("A node recorder.\n");
}
