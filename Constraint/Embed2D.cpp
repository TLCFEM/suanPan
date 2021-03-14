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
	: Constraint(T, S, 0, {NT}, {ET}, 2) {}

int Embed2D::initialize(const shared_ptr<DomainBase>& D) {
	auto& t_node = D->get<Node>(node_encoding(0));
	auto& t_element = D->get<Element>(dof_reference(0));

	if(nullptr == t_node || nullptr == t_element || !t_node->is_active() || !t_element->is_active() || 4 != t_element->get_node_number()) {
		D->disable_constraint(get_tag());
		return SUANPAN_FAIL;
	}

	const auto t_coor = get_coordinate(t_element.get(), 2);

	const vec n_coor = t_node->get_coordinate().head(2);

	vec t_para = zeros(2);

	rowvec n;

	unsigned counter = 0;
	while(true) {
		if(max_iteration == ++counter) {
			D->disable_constraint(get_tag());
			return SUANPAN_FAIL;
		}

		const vec incre = solve((shape::quad(t_para, 1, 4) * t_coor).t(), n_coor - ((n = shape::quad(t_para, 0, 4)) * t_coor).t());
		if(norm(incre) < 1E-14) break;
		t_para += incre;
	}

	auto& n_dof = t_node->get_reordered_dof();
	auto& e_dof = t_element->get_dof_encoding();

	auxiliary_stiffness.zeros(D->get_factory()->get_size(), num_size);

	for(auto K = 0u; K < num_size; ++K) {
		auxiliary_stiffness(n_dof(K), K) = -1.;
		for(uword I = 0, J = K; I < n.n_elem; ++I, J += num_size) auxiliary_stiffness(e_dof(J), K) = n(I);
	}

	return Constraint::initialize(D);
}

int Embed2D::process(const shared_ptr<DomainBase>& D) {
	auxiliary_resistance = auxiliary_stiffness.t() * D->get_factory()->get_trial_displacement();

	return SUANPAN_SUCCESS;
}
