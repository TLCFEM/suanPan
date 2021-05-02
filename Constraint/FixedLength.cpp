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

#include "FixedLength.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

FixedLength::FixedLength(const unsigned T, const unsigned S, const unsigned A, const unsigned D, uvec&& N)
	: Constraint(T, S, A, std::forward<uvec>(N), 2 == D ? uvec{1, 2} : uvec{1, 2, 3}, 1)
	, inequal(false)
	, min_gap(0.) { set_connected(true); }

FixedLength::FixedLength(const unsigned T, const unsigned S, const unsigned A, const unsigned D, const double M, uvec&& N)
	: Constraint(T, S, A, std::forward<uvec>(N), 2 == D ? uvec{1, 2} : uvec{1, 2, 3}, 1)
	, inequal(true)
	, min_gap(M * M) { set_connected(true); }

int FixedLength::initialize(const shared_ptr<DomainBase>& D) {
	dof_encoding = get_nodal_active_dof(D);

	if(dof_encoding.n_elem != node_encoding.n_elem * dof_reference.n_elem) {
		D->disable_constraint(get_tag());
		return SUANPAN_SUCCESS;
	}

	coor = resize(D->get<Node>(node_encoding(1))->get_coordinate(), dof_reference.n_elem, 1) - resize(D->get<Node>(node_encoding(0))->get_coordinate(), dof_reference.n_elem, 1);

	current_resistance = trial_resistance.zeros(num_size);

	return Constraint::initialize(D);
}

int FixedLength::process(const shared_ptr<DomainBase>& D) {
	auto& W = D->get_factory();

	const auto& n_dof = dof_reference.n_elem;

	const uvec dof_i = dof_encoding.head(n_dof);
	const uvec dof_j = dof_encoding.tail(n_dof);

	const vec t_disp = W->get_trial_displacement()(dof_j) - W->get_trial_displacement()(dof_i);

	if(inequal) {
		if(accu(square(coor + t_disp)) > min_gap + datum::eps) {
			set_multiplier_size(0);
			return SUANPAN_SUCCESS;
		}
		set_multiplier_size(1);
		auxiliary_load = min_gap - dot(coor, coor);
	}

	auxiliary_stiffness.zeros(W->get_size(), num_size);
	auxiliary_resistance = 0.;
	for(auto I = 0llu; I < n_dof; ++I) {
		auxiliary_stiffness(dof_i(I)) = -(auxiliary_stiffness(dof_j(I)) = 2. * (coor(I) + t_disp(I)));
		auxiliary_resistance += t_disp(I) * (2. * coor(I) + t_disp(I));
	}

	stiffness.zeros(dof_encoding.n_elem, dof_encoding.n_elem);
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

MinimumGap::MinimumGap(const unsigned T, const unsigned S, const unsigned A, const unsigned D, const double M, uvec&& N)
	: FixedLength(T, S, A, D, M, std::forward<uvec>(N)) {}
