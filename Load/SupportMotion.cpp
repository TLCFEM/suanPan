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

#include "SupportMotion.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Load/Amplitude/Amplitude.h>
#include <Solver/Integrator/Integrator.h>
#include <Step/Step.h>

SupportMotion::SupportMotion(const unsigned T, const unsigned ST, const double L, uvec&& N, const unsigned D, const unsigned AT)
	: Load(T, ST, AT, std::forward<uvec>(N), uvec{D}, L) { enable_displacement_control(); }

SupportMotion::SupportMotion(const unsigned T, const unsigned ST, const double L, uvec&& N, uvec&& D, const unsigned AT)
	: Load(T, ST, AT, std::forward<uvec>(N), std::forward<uvec>(D), L) { enable_displacement_control(); }

int SupportMotion::initialize(const shared_ptr<DomainBase>& D) {
	if(initialized) return SUANPAN_SUCCESS;

	if(SUANPAN_SUCCESS != Load::initialize(D)) return SUANPAN_FAIL;

	set_end_step(start_step + 1);

	const auto& W = D->get_factory();

	vector<uword> r_dof, e_dof;

	for(auto I : W->get_reference_dof()) r_dof.emplace_back(I);

	for(auto I : node_encoding)
		if(auto& t_node = D->get<Node>(I); t_node != nullptr && t_node->is_active()) {
			auto& t_dof = t_node->get_reordered_dof();
			for(const auto J : dof_reference)
				if(J < t_dof.n_elem) {
					e_dof.emplace_back(t_dof(J));
					if(find(r_dof.begin(), r_dof.end(), e_dof.back()) == r_dof.end()) r_dof.emplace_back(e_dof.back());
					else suanpan_warning("more than one displacement loads are applied on node %llu DoF %llu.\n", I, J);
				}
		}

	encoding = e_dof;

	W->set_reference_dof(uvec(r_dof));
	W->set_reference_size(static_cast<unsigned>(r_dof.size()));

	return SUANPAN_SUCCESS;
}

int SupportDisplacement::process(const shared_ptr<DomainBase>& D) {
	const auto& W = D->get_factory();

	trial_settlement.zeros(W->get_size());

	trial_settlement(encoding).fill(pattern * magnitude->get_amplitude(W->get_trial_time()));

	return SUANPAN_SUCCESS;
}

int SupportVelocity::process(const shared_ptr<DomainBase>& D) {
	const auto& W = D->get_factory();
	const auto& G = D->get_current_step()->get_integrator();

	trial_settlement.zeros(W->get_size());

	trial_settlement(encoding) = G->from_total_velocity(pattern * magnitude->get_amplitude(W->get_trial_time()), encoding);

	return SUANPAN_SUCCESS;
}

int SupportAcceleration::process(const shared_ptr<DomainBase>& D) {
	const auto& W = D->get_factory();
	const auto& G = D->get_current_step()->get_integrator();

	trial_settlement.zeros(W->get_size());

	trial_settlement(encoding) = G->from_total_acceleration(pattern * magnitude->get_amplitude(W->get_trial_time()), encoding);

	return SUANPAN_SUCCESS;
}
