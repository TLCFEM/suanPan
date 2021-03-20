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

#include "FixedLength.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

FixedLength::FixedLength(const unsigned T, const unsigned S, const unsigned A, const unsigned D, uvec&& N)
	: Constraint(T, S, A, std::forward<uvec>(N), {D}, 1) {}

int FixedLength::initialize(const shared_ptr<DomainBase>& D) {
	for(uword I = 0; I < node_encoding.n_elem; ++I) {
		auto& t_node = D->get<Node>(node_encoding(I));
		if(nullptr == t_node || !t_node->is_active() || t_node->get_reordered_dof().n_elem < dof_reference(0)) {
			D->disable_constraint(get_tag());
			return SUANPAN_SUCCESS;
		}
	}

	set_connected(true);

	node_i = D->get<Node>(node_encoding(0));
	node_j = D->get<Node>(node_encoding(1));

	const auto& n_dof = dof_reference(0);

	dof_encoding = join_cols(node_i.lock()->get_reordered_dof().head(n_dof), node_j.lock()->get_reordered_dof().head(n_dof));

	current_resistance = trial_resistance.zeros(num_size);

	return Constraint::initialize(D);
}

int FixedLength::process(const shared_ptr<DomainBase>& D) {
	const auto& node_ptr_i = node_i.lock();
	const auto& node_ptr_j = node_j.lock();

	const auto& n_dof = dof_reference(0);

	vec coor = resize(node_ptr_j->get_coordinate(), n_dof, 1) - resize(node_ptr_i->get_coordinate(), n_dof, 1);
	vec t_disp = node_ptr_j->get_trial_displacement().head(n_dof) - node_ptr_i->get_trial_displacement().head(n_dof);
	uvec dof_i = node_ptr_i->get_reordered_dof().head(n_dof);
	uvec dof_j = node_ptr_j->get_reordered_dof().head(n_dof);

	auto& W = D->get_factory();

	auxiliary_stiffness.zeros(W->get_size(), num_size);
	auxiliary_resistance.zeros();
	for(auto I = 0llu; I < n_dof; ++I) {
		auxiliary_stiffness(dof_i(I)) = -(auxiliary_stiffness(dof_j(I)) = 2. * (coor(I) + t_disp(I)));
		auxiliary_resistance += t_disp(I) * (2. * coor(I) + t_disp(I));
	}

	stiffness.zeros(2 * n_dof, 2 * n_dof);
	const auto t_factor = 2. * trial_lambda(0);
	for(auto I = 0llu; I < n_dof; ++I) stiffness(I + n_dof, I) = stiffness(I, I + n_dof) = -(stiffness(I, I) = stiffness(I + n_dof, I + n_dof) = t_factor);

	trial_resistance = auxiliary_stiffness * trial_lambda;

	return SUANPAN_SUCCESS;
}

void FixedLength::update_status(const vec& i_lambda) { trial_lambda += i_lambda; }

void FixedLength::commit_status() {
	current_lambda = trial_lambda;
	current_resistance = trial_resistance;
}

void FixedLength::clear_status() {
	current_lambda = trial_lambda.zeros();
	current_resistance = trial_resistance.zeros();
}

void FixedLength::reset_status() {
	trial_lambda = current_lambda;
	trial_resistance = current_resistance;
}
