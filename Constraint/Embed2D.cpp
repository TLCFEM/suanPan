////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2021 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "Embed2D.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Element/Element.h>
#include <Toolbox/shapeFunction.h>

Embed2D::Embed2D(const unsigned T, const unsigned S, const unsigned ET, const unsigned NT)
	: Constraint(T, S, 0, {ET, NT}, {}) {}

int Embed2D::initialize(const shared_ptr<DomainBase>& D) {
	if(!D->find<Node>(nodes(1)) || !D->find<Element>(nodes(0))) {
		D->disable_constraint(get_tag());
		return SUANPAN_SUCCESS;
	}

	auto& t_node = D->get<Node>(nodes(1));
	auto& t_element = D->get<Element>(nodes(0));

	if(!t_node->is_active() || !t_element->is_active() || 4 != t_element->get_node_number()) {
		D->disable_constraint(get_tag());
		return SUANPAN_SUCCESS;
	}

	const auto t_coor = get_coordinate(t_element.get(), 2);

	const vec n_coor = t_node->get_coordinate().head(2);

	vec t_para = zeros(2);

	auto& n = access::rw(weight);

	unsigned counter = 0;
	while(true) {
		if(max_iteration == ++counter) {
			D->disable_constraint(get_tag());
			return SUANPAN_SUCCESS;
		}

		const vec incre = solve((shape::quad(t_para, 1, 4) * t_coor).t(), n_coor - ((n = shape::quad(t_para, 0, 4)) * t_coor).t());
		if(norm(incre) < 1E-14) break;
		t_para += incre;
	}

	return Constraint::initialize(D);
}

int Embed2D::process(const shared_ptr<DomainBase>& D) {
	auto& W = D->get_factory();

	auto last_pos = W->get_mpc();

	auto& n_dof = D->get<Node>(nodes(1))->get_reordered_dof();
	auto& e_dof = D->get<Element>(nodes(0))->get_dof_encoding();

	sp_vec slice_x(W->get_size());
	slice_x(n_dof(0)) = -1.;
	for(uword I = 0, J = 0; I < weight.n_elem; ++I, J += 2) slice_x(e_dof(J)) = weight(I);

	W->incre_mpc();
	get_auxiliary_stiffness(W).col(last_pos++) = slice_x;
	// get_auxiliary_load(W).back() = 0.;

	sp_vec slice_y(W->get_size());
	slice_y(n_dof(1)) = -1.;
	for(uword I = 0, J = 1; I < weight.n_elem; ++I, J += 2) slice_y(e_dof(J)) = weight(I);

	W->incre_mpc();
	get_auxiliary_stiffness(W).col(last_pos) = slice_y;
	// get_auxiliary_load(W).back() = 0.;

	return SUANPAN_SUCCESS;
}
