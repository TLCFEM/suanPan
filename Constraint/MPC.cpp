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

#include "MPC.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

MPC::MPC(const unsigned T, const unsigned S, const unsigned A, uvec&& N, uvec&& D, vec&& W, const double L)
	: Constraint(T, S, A, std::forward<uvec>(N), std::forward<uvec>(D))
	, weight_pool(std::forward<vec>(W))
	, psudo_load(L) {}

int MPC::initialize(const shared_ptr<DomainBase>& D) {
	auto flag = true;
	for(uword I = 0; I < nodes.n_elem; ++I) {
		if(!D->find<Node>(nodes(I))) {
			flag = false;
			break;
		}
		if(auto& t_node = D->get<Node>(nodes(I)); !t_node->is_active() || t_node->get_reordered_dof().n_elem <= dofs(I)) {
			flag = false;
			break;
		}
	}

	if(!flag) {
		D->disable_constraint(get_tag());
		return SUANPAN_FAIL;
	}

	return Constraint::initialize(D);
}

int MPC::process(const shared_ptr<DomainBase>& D) {
	auto& W = D->get_factory();

	const auto last_pos = W->get_mpc();

	auto flag = true;
	auto resistance = 0.;
	sp_vec slice(W->get_size());
	for(uword I = 0; I < nodes.n_elem; ++I) {
		auto& t_node = D->get<Node>(nodes(I));
		if(nullptr == t_node || !t_node->is_active()) {
			flag = false;
			break;
		}
		auto& t_dof = t_node->get_reordered_dof();
		auto& t_disp = t_node->get_trial_displacement();
		if(t_dof.n_elem <= dofs(I)) {
			flag = false;
			break;
		}
		slice(t_dof(dofs(I) - 1)) = weight_pool(I);
		resistance += weight_pool(I) * t_disp(dofs(I) - 1);
	}

	if(flag) {
		W->incre_mpc();
		get_auxiliary_stiffness(W).col(last_pos) = slice;
		get_auxiliary_load(W).back() = psudo_load * magnitude->get_amplitude(W->get_trial_time());
		get_trial_auxiliary_resistance(W).back() = resistance;
	}
	else suanpan_debug("some node or DoF is not active, the MPC %u is not applied.\n", get_tag());

	return SUANPAN_SUCCESS;
}
