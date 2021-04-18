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

#include "MPC.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Load/Amplitude/Amplitude.h>

MPC::MPC(const unsigned T, const unsigned S, const unsigned A, uvec&& N, uvec&& D, vec&& W, const double L)
	: Constraint(T, S, A, std::forward<uvec>(N), std::forward<uvec>(D), 1)
	, weight_pool(std::forward<vec>(W))
	, psudo_load(L) {}

int MPC::initialize(const shared_ptr<DomainBase>& D) {
	auto& W = D->get_factory();

	auxiliary_stiffness.zeros(W->get_size(), num_size);

	for(uword I = 0; I < node_encoding.n_elem; ++I) {
		auto& t_node = D->get<Node>(node_encoding(I));
		if(nullptr == t_node || !t_node->is_active() || t_node->get_reordered_dof().n_elem <= dof_reference(I)) {
			auxiliary_stiffness.reset();
			D->disable_constraint(get_tag());
			return SUANPAN_SUCCESS;
		}
		auto& t_dof = t_node->get_reordered_dof();
		auxiliary_stiffness(t_dof(dof_reference(I) - 1)) = weight_pool(I);
	}

	return Constraint::initialize(D);
}

int MPC::process(const shared_ptr<DomainBase>& D) {
	auto& W = D->get_factory();

	auxiliary_load = psudo_load * magnitude->get_amplitude(W->get_trial_time());

	auxiliary_resistance = auxiliary_stiffness.t() * W->get_trial_displacement();

	return SUANPAN_SUCCESS;
}
