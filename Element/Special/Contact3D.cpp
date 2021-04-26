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

#include "Contact3D.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Group.h>
#include <Domain/Node.h>

void Contact3D::update_position() {
	for(auto& [node, position, ne, nf] : master) {
		for(auto I = 0; I < 3; ++I) {
			position[I].zeros(3);
			auto& t_coor = node[I].lock()->get_coordinate();
			for(uword K = 0; K < std::min(3llu, t_coor.n_elem); ++K) position[I](K) += t_coor(K);
			auto& t_disp = node[I].lock()->get_trial_displacement();
			for(uword K = 0; K < std::min(3llu, t_disp.n_elem); ++K) position[I](K) += t_disp(K);
		}

		const vec edge_i = position[1] - position[0];
		const vec edge_j = position[2] - position[1];
		const vec edge_k = position[0] - position[2];

		nf = normalise(cross(edge_i, edge_j));

		ne[0] = normalise(cross(edge_i, nf));
		ne[1] = normalise(cross(edge_j, nf));
		ne[2] = normalise(cross(edge_k, nf));
	}

	for(auto& [node, position] : slave) {
		position.zeros(3);
		auto& t_coor = node.lock()->get_coordinate();
		for(uword K = 0; K < std::min(3llu, t_coor.n_elem); ++K) position(K) += t_coor(K);
		auto& t_disp = node.lock()->get_trial_displacement();
		for(uword K = 0; K < std::min(3llu, t_disp.n_elem); ++K) position(K) += t_disp(K);
	}
}

Contact3D::Contact3D(const unsigned T, const unsigned M, const unsigned S, const double P)
	: Element(T, c_dof, {M, S})
	, master_tag(M)
	, slave_tag(S)
	, alpha(P) {}

void Contact3D::initialize(const shared_ptr<DomainBase>& D) {
	if(!D->find<Group>(master_tag) || !D->find<Group>(slave_tag)) {
		suanpan_error("Contact3D: cannot find the given groups for element %u.\n", get_tag());
		D->disable_element(get_tag());
		return;
	}

	const auto& m_pool = D->get<Group>(master_tag)->get_pool();

	if(0 != m_pool.n_elem % 3) {
		suanpan_error("Contact3D: master group has wrong number of nodes, now disable element %u and return.", get_tag());
		D->disable_element(get_tag());
		return;
	}

	master.reserve(m_pool.n_elem / 3);
	for(auto I = 0llu, J = 1llu, K = 2llu; K < m_pool.n_elem; I += 3llu, J += 3llu, K += 3llu) {
		if(!D->find<Node>(m_pool(I)) || !D->find<Node>(m_pool(J)) || !D->find<Node>(m_pool(K))) {
			suanpan_error("Contact3D: invalid Node %llu detected, now disable element %u and return.", I, get_tag());
			D->disable_element(get_tag());
			return;
		}
		master.emplace_back(MasterFacet{{D->get<Node>(m_pool(I)), D->get<Node>(m_pool(J)), D->get<Node>(m_pool(K))}, {}, {}, {}});
	}
	master.shrink_to_fit();

	const auto& s_pool = D->get<Group>(slave_tag)->get_pool();
	slave.reserve(s_pool.n_elem);
	for(auto& I : s_pool) {
		if(!D->find<Node>(I)) {
			suanpan_error("Contact3D: invalid Node %llu detected, now disable element %u and return.", I, get_tag());
			D->disable_element(get_tag());
			return;
		}
		slave.emplace_back(SlaveNode{D->get<Node>(I), {}});
	}
	slave.shrink_to_fit();
}

int Contact3D::update_status() {
	update_position();

	trial_stiffness.zeros(get_total_number(), get_total_number());
	trial_resistance.zeros(get_total_number());

	return SUANPAN_SUCCESS;
}

int Contact3D::clear_status() {
	current_stiffness.reset();
	trial_stiffness.reset();
	current_resistance.reset();
	trial_resistance.reset();
	return SUANPAN_SUCCESS;
}

int Contact3D::commit_status() {
	current_stiffness = trial_stiffness;
	current_resistance = trial_resistance;
	return SUANPAN_SUCCESS;
}

int Contact3D::reset_status() {
	trial_stiffness = current_stiffness;
	trial_resistance = current_resistance;
	return SUANPAN_SUCCESS;
}
