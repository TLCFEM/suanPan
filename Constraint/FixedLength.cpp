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
	auto flag = true;
	for(uword I = 0; I < nodes.n_elem; ++I) {
		if(!D->find<Node>(nodes(I))) {
			flag = false;
			break;
		}
		if(auto& t_node = D->get<Node>(nodes(I)); !t_node->is_active() || t_node->get_reordered_dof().n_elem <= dofs(0)) {
			flag = false;
			break;
		}
	}

	if(flag) return Constraint::initialize(D);

	D->disable_constraint(get_tag());
	return SUANPAN_FAIL;
}

int FixedLength::process(const shared_ptr<DomainBase>& D) {
	auto& W = D->get_factory();

	auto& node_i = D->get<Node>(nodes(0));
	auto& node_j = D->get<Node>(nodes(1));

	vec coor_i = resize(node_i->get_coordinate(), dofs(0), 1);
	vec coor_j = resize(node_j->get_coordinate(), dofs(0), 1);
	vec t_disp_i = node_i->get_trial_displacement().head(dofs(0));
	vec t_disp_j = node_j->get_trial_displacement().head(dofs(0));
	vec c_disp_i = node_i->get_current_displacement().head(dofs(0));
	vec c_disp_j = node_j->get_current_displacement().head(dofs(0));
	uvec dof_i = node_i->get_reordered_dof().head(dofs(0));
	uvec dof_j = node_j->get_reordered_dof().head(dofs(0));

	const auto last_pos = W->get_mpc();

	W->incre_mpc();

	auto c_resistance = 0., t_resistance = 0.;
	sp_vec slice(W->get_size());

	for(auto I = 0llu; I < dofs(0); ++I) {
		slice(dof_i(I)) = -(slice(dof_j(I)) = 2. * (coor_j(I) - coor_i(I) + t_disp_j(I) - t_disp_i(I)));
		c_resistance += (c_disp_j(I) - c_disp_i(I)) * (2. * (coor_j(I) - coor_i(I)) + c_disp_j(I) - c_disp_i(I));
		t_resistance += (t_disp_j(I) - t_disp_i(I)) * (2. * (coor_j(I) - coor_i(I)) + t_disp_j(I) - t_disp_i(I));
	}

	get_auxiliary_stiffness(W).col(last_pos) = slice;
	get_current_auxiliary_resistance(W).back() = c_resistance;
	get_trial_auxiliary_resistance(W).back() = t_resistance;
	get_auxiliary_encoding(W).back() = get_tag();

	auto& t_matrix = W->get_stiffness();
	auto t_factor = 2. * trial_lambda(0);
	if(StorageScheme::SPARSE == D->get_factory()->get_storage_scheme())
		for(auto I = 0llu; I < dofs(0); ++I) {
			t_matrix->at(dof_i(I), dof_i(I)) = t_factor;
			t_matrix->at(dof_j(I), dof_j(I)) = t_factor;
			t_matrix->at(dof_i(I), dof_j(I)) = -t_factor;
			t_matrix->at(dof_j(I), dof_i(I)) = -t_factor;
		}
	else
		for(auto I = 0llu; I < dofs(0); ++I) {
			t_matrix->at(dof_i(I), dof_i(I)) += t_factor;
			t_matrix->at(dof_j(I), dof_j(I)) += t_factor;
			t_matrix->at(dof_i(I), dof_j(I)) -= t_factor;
			t_matrix->at(dof_j(I), dof_i(I)) -= t_factor;
		}

	W->update_trial_load(W->get_trial_load() - trial_lambda(0) * slice);

	return SUANPAN_SUCCESS;
}
