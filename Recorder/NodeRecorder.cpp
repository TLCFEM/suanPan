/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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
	for(const auto I : get_object_tag())
		if(!D->find<Node>(I)) {
			D->disable_recorder(get_tag());
			return;
		}
}

void NodeRecorder::record(const shared_ptr<DomainBase>& D) {
	if(1 != interval && counter++ != interval) return;

	counter = 1;

	auto& obj_tag = get_object_tag();

	auto insert_damping_force = [&](const uword J) {
		for(unsigned I = 0; I < obj_tag.n_elem; ++I) {
			const auto& t_dof = D->get<Node>(obj_tag(I))->get_reordered_dof();
			const auto& t_force = D->get_factory()->get_current_damping_force();
			insert({{(t_dof.n_elem > J && t_force.n_elem > t_dof(J) ? t_force(t_dof(J)) : 0.)}}, I);
		}
	};
	auto insert_inertial_force = [&](const uword J) {
		for(unsigned I = 0; I < obj_tag.n_elem; ++I) {
			const auto& t_dof = D->get<Node>(obj_tag(I))->get_reordered_dof();
			const auto& t_force = D->get_factory()->get_current_inertial_force();
			insert({{(t_dof.n_elem > J && t_force.n_elem > t_dof(J) ? t_force(t_dof(J)) : 0.)}}, I);
		}
	};

	if(OutputType::DF == get_variable_type()) for(unsigned I = 0; I < obj_tag.n_elem; ++I) insert({D->get_factory()->get_current_damping_force()(D->get<Node>(obj_tag(I))->get_reordered_dof())}, I);
	else if(OutputType::DF1 == get_variable_type()) insert_damping_force(0);
	else if(OutputType::DF2 == get_variable_type()) insert_damping_force(1);
	else if(OutputType::DF3 == get_variable_type()) insert_damping_force(2);
	else if(OutputType::DF4 == get_variable_type() || OutputType::DM1 == get_variable_type()) insert_damping_force(3);
	else if(OutputType::DF5 == get_variable_type() || OutputType::DM2 == get_variable_type()) insert_damping_force(4);
	else if(OutputType::DF6 == get_variable_type() || OutputType::DM3 == get_variable_type()) insert_damping_force(5);
	else if(OutputType::IF == get_variable_type()) for(unsigned I = 0; I < obj_tag.n_elem; ++I) insert({D->get_factory()->get_current_inertial_force()(D->get<Node>(obj_tag(I))->get_reordered_dof())}, I);
	else if(OutputType::IF1 == get_variable_type()) insert_inertial_force(0);
	else if(OutputType::IF2 == get_variable_type()) insert_inertial_force(1);
	else if(OutputType::IF3 == get_variable_type()) insert_inertial_force(2);
	else if(OutputType::IF4 == get_variable_type() || OutputType::IM1 == get_variable_type()) insert_inertial_force(3);
	else if(OutputType::IF5 == get_variable_type() || OutputType::IM2 == get_variable_type()) insert_inertial_force(4);
	else if(OutputType::IF6 == get_variable_type() || OutputType::IM3 == get_variable_type()) insert_inertial_force(5);
	else for(unsigned I = 0; I < obj_tag.n_elem; ++I) insert(D->get<Node>(obj_tag(I))->record(get_variable_type()), I);

	if(if_record_time()) insert(D->get_factory()->get_current_time());
}

void NodeRecorder::print() { suanpan_info("A Node Recorder.\n"); }
