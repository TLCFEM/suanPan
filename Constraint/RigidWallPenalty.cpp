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

#include "RigidWallPenalty.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

RigidWallPenalty::RigidWallPenalty(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& N, const double F)
	: Constraint(T, S, A, {}, {}, 0)
	, alpha(F)
	, origin(std::forward<vec>(O))
	, norm(normalise(N)) {}

RigidWallPenalty::RigidWallPenalty(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& E1, vec&& E2, const double F)
	: Constraint(T, S, A, {}, {}, 0)
	, alpha(F)
	, edge_a(E1 - O)
	, edge_b(E2 - O)
	, origin(std::forward<vec>(O))
	, norm(normalise(cross(E1, E2))) {}

int RigidWallPenalty::initialize(const shared_ptr<DomainBase>& D) {
	if(origin.n_elem != norm.n_elem || !edge_a.empty() && edge_a.n_elem != norm.n_elem || !edge_b.empty() && edge_b.n_elem != norm.n_elem) {
		D->disable_constraint(get_tag());
		return SUANPAN_SUCCESS;
	}

	return Constraint::initialize(D);
}

int RigidWallPenalty::process(const shared_ptr<DomainBase>& D) {
	auto& W = D->get_factory();
	auto& t_stiff = W->get_stiffness();
	auto& t_load = get_trial_load(W);

	const auto factor = alpha * pow(W->get_incre_time(), -2.);

	// penalty method
	if(StorageScheme::SPARSE == W->get_storage_scheme())
		for(const auto& I : D->get_node_pool()) {
			auto& t_dof = I->get_reordered_dof();
			auto& t_coor = I->get_coordinate();
			auto& t_disp = I->get_trial_displacement();
			const auto size = std::min(norm.n_elem, std::min(t_coor.n_elem, t_disp.n_elem));
			vec new_pos(norm.n_elem, fill::zeros);
			for(auto J = 0llu; J < size; ++J) new_pos(J) = t_coor(J) + t_disp(J);
			new_pos -= origin;
			if(!edge_a.empty() && dot(new_pos, edge_a) > arma::norm(edge_a)) continue;
			if(!edge_b.empty() && dot(new_pos, edge_b) > arma::norm(edge_b)) continue;
			const auto pen = dot(new_pos, norm);
			if(pen >= datum::eps) continue;
			for(auto J = 0llu; J < size; ++J) {
				for(auto K = 0llu; K < size; ++K) t_stiff->at(t_dof(J), t_dof(K)) = factor * norm(J) * norm(K);
				t_load(t_dof(J)) -= factor * pen * norm(J);
			}
		}
	else
		for(const auto& I : D->get_node_pool()) {
			auto& t_dof = I->get_reordered_dof();
			auto& t_coor = I->get_coordinate();
			auto& t_disp = I->get_trial_displacement();
			const auto size = std::min(norm.n_elem, std::min(t_coor.n_elem, t_disp.n_elem));
			vec new_pos(norm.n_elem, fill::zeros);
			for(auto J = 0llu; J < size; ++J) new_pos(J) = t_coor(J) + t_disp(J);
			new_pos -= origin;
			if(!edge_a.empty() && dot(new_pos, edge_a) > arma::norm(edge_a)) continue;
			if(!edge_b.empty() && dot(new_pos, edge_b) > arma::norm(edge_b)) continue;
			const auto pen = dot(new_pos, norm);
			if(pen >= datum::eps) continue;
			for(auto J = 0llu; J < size; ++J) {
				for(auto K = 0llu; K < size; ++K) t_stiff->at(t_dof(J), t_dof(K)) += factor * norm(J) * norm(K);
				t_load(t_dof(J)) -= factor * pen * norm(J);
			}
		}

	return SUANPAN_SUCCESS;
}
