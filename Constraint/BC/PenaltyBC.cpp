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

#include "PenaltyBC.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Toolbox/utility.h>

/**
 * \brief default constructor.
 * \param T `unique_tag`
 * \param S `start_step`
 * \param N `nodes`
 * \param D `dofs`
 */
PenaltyBC::PenaltyBC(const unsigned T, const unsigned S, uvec&& N, const unsigned D)
	: Constraint(T, S, 0, std::forward<uvec>(N), {D}, 0) {}

/**
 * \brief the constructor uses DoF vector.
 * \param T `unique_tag`
 * \param S `start_step`
 * \param N `nodes`
 * \param D `dofs`
 */
PenaltyBC::PenaltyBC(const unsigned T, const unsigned S, uvec&& N, uvec&& D)
	: Constraint(T, S, 0, std::forward<uvec>(N), std::forward<uvec>(D), 0) {}

/**
 * \brief the constructor uses predefined TYPE: "XSYMM", "YSYMM",
 * "ZSYMM", "ENCASTRE",
 * "PINNED".
 * \param T `unique_tag`
 * \param S `start_step`
 * \param N `nodes`
 * \param TP PenaltyBC TYPE
 */
PenaltyBC::PenaltyBC(const unsigned T, const unsigned S, uvec&& N, const char* TP)
	: Constraint(T, S, 0, std::forward<uvec>(N), {}, 0) {
	if(is_equal(TP[0], 'X')) dofs = uvec{1, 5, 6};
	else if(is_equal(TP[0], 'Y')) dofs = uvec{2, 4, 6};
	else if(is_equal(TP[0], 'Z')) dofs = uvec{3, 4, 5};
	else if(is_equal(TP[0], 'E')) dofs = uvec{1, 2, 3, 4, 5, 6};
	else if(is_equal(TP[0], 'P')) dofs = uvec{1, 2, 3};
	else if(is_equal(TP[0], '1')) dofs = uvec{1};
	else if(is_equal(TP[0], '2')) dofs = uvec{2};
	else if(is_equal(TP[0], '3')) dofs = uvec{3};
	else if(is_equal(TP[0], '4')) dofs = uvec{4};
	else if(is_equal(TP[0], '5')) dofs = uvec{5};
	else if(is_equal(TP[0], '6')) dofs = uvec{6};
}

/**
 * \brief default destructor.
 */
PenaltyBC::~PenaltyBC() = default;

/**
 * \brief method to apply the PenaltyBC to the system.
 * \param D `Domain`
 * \return 0
 */
int PenaltyBC::process(const shared_ptr<DomainBase>& D) {
	auto& t_matrix = D->get_factory()->get_stiffness();

	if(StorageScheme::SPARSE == D->get_factory()->get_storage_scheme()) {
		const auto max_term = std::min(multiplier * t_matrix->max(), 1E13);

		for(const auto& I : nodes)
			if(D->find<Node>(I))
				if(auto& t_node = D->get<Node>(I); t_node->is_active()) {
					auto& t_dof = t_node->get_reordered_dof();
					for(const auto& J : dofs) if(J <= t_dof.n_elem) if(auto& t_idx = t_dof(J - 1); D->insert_restrained_dof(t_idx)) t_matrix->at(t_idx, t_idx) = max_term;
				}
	}
	else {
		auto& t_set_b = D->get_constrained_dof();

		for(const auto& I : nodes)
			if(D->find<Node>(I))
				if(auto& t_node = D->get<Node>(I); t_node->is_active()) {
					auto& t_dof = t_node->get_reordered_dof();
					for(const auto& J : dofs)
						if(J <= t_dof.n_elem)
							if(auto& t_idx = t_dof(J - 1); D->insert_restrained_dof(t_idx)) {
								if(0. != t_matrix->operator()(t_idx, t_idx)) t_matrix->at(t_idx, t_idx) = std::min(multiplier * t_matrix->operator()(t_idx, t_idx), 1E13);
								else {
									auto& t_set = D->get_restrained_dof();
									t_matrix->at(t_idx, t_idx) = t_set.size() == 1 ? t_set_b.empty() ? std::min(multiplier * t_matrix->max(), 1E13) : t_matrix->operator()(*t_set_b.cbegin(), *t_set_b.cbegin()) : *t_set.cbegin() == t_idx ? t_matrix->operator()(*++t_set.cbegin(), *++t_set.cbegin()) : t_matrix->operator()(*t_set.cbegin(), *t_set.cbegin());
								}
							}
				}
	}

	return SUANPAN_SUCCESS;
}
