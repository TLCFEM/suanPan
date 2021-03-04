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

#include "RigidWallMultiplier.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

int RigidWallMultiplier::initialize(const shared_ptr<DomainBase>& D) {
	// for dynamic analysis, velocity and acceleration shall be changed
	// multiplier method does not give proper result
	if(AnalysisType::DYNAMICS == D->get_factory()->get_analysis_type()) {
		suanpan_warning("multiplier rigid wall constraint cannot be applied in dynamic analysis.\n");
		access::rw(use_penalty) = true;
	}

	return RigidWallPenalty::initialize(D);
}

int RigidWallMultiplier::process(const shared_ptr<DomainBase>& D) {
	if(use_penalty) return RigidWallPenalty::process(D);

	auto& W = D->get_factory();

	// multiplier method
	for(const auto& I : D->get_node_pool()) {
		auto& t_dof = I->get_reordered_dof();
		auto& t_coor = I->get_coordinate();
		auto& t_disp = I->get_trial_displacement();
		const auto size = std::min(std::min(t_coor.n_elem, t_disp.n_elem), norm.n_elem);
		vec new_pos(norm.n_elem, fill::zeros);
		for(auto J = 0llu; J < size; ++J) new_pos(J) = t_coor(J) + t_disp(J);
		new_pos -= origin;
		if(!edge_a.empty() && dot(new_pos, edge_a) > arma::norm(edge_a)) continue;
		if(!edge_b.empty() && dot(new_pos, edge_b) > arma::norm(edge_b)) continue;
		if(const auto pen = dot(new_pos, norm); pen < datum::eps)
			for(auto J = 0llu; J < size; ++J) {
				W->incre_mpc();
				get_auxiliary_stiffness(W)(t_dof(J), W->get_mpc() - 1) = -1.;
				get_trial_auxiliary_resistance(W).back() = -pen * norm(J);
			}
	}

	return SUANPAN_SUCCESS;
}
