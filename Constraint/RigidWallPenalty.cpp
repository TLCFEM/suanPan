/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "RigidWallPenalty.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

RigidWallPenalty::RigidWallPenalty(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& N, const double F)
    : Constraint(T, S, A, {}, {}, 0)
    , alpha(F)
    , origin(std::forward<vec>(O))
    , outer_norm(normalise(N)) {}

RigidWallPenalty::RigidWallPenalty(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& E1, vec&& E2, const double F)
    : Constraint(T, S, A, {}, {}, 0)
    , alpha(F)
    , edge_a(E1 - O)
    , edge_b(E2 - O)
    , origin(std::forward<vec>(O))
    , outer_norm(normalise(cross(E1, E2)))
    , length_a(norm(edge_a))
    , length_b(norm(edge_b)) {}

int RigidWallPenalty::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    const auto factor = alpha * pow(W->get_incre_time(), -2.);

    stiffness.reset();
    resistance.zeros(W->get_size());
    vector<uword> pool;

    auto counter = 0llu;
    for(const auto& I : D->get_node_pool()) {
        auto& t_coor = I->get_coordinate();
        auto& t_disp = I->get_trial_displacement();
        const auto t_size = std::min(outer_norm.n_elem, std::min(t_coor.n_elem, t_disp.n_elem));
        vec t_pos = -origin;
        for(auto J = 0llu; J < t_size; ++J) t_pos(J) += t_coor(J) + t_disp(J);
        if(!edge_a.empty() && dot(t_pos, edge_a) > length_a || !edge_b.empty() && dot(t_pos, edge_b) > length_b) continue;
        const auto t_pen = dot(t_pos, outer_norm);
        if(t_pen > datum::eps) continue;
        const auto next_counter = counter + t_size;
        stiffness.resize(next_counter, next_counter);
        stiffness.submat(counter, counter, size(t_size, t_size)) = factor * outer_norm.head(t_size) * outer_norm.head(t_size).t();
        auto& t_dof = I->get_reordered_dof();
        for(auto J = 0llu; J < t_size; ++J) {
            pool.emplace_back(t_dof(J));
            resistance(t_dof(J)) += factor * t_pen * outer_norm(J);
        }
        counter = next_counter;
    }

    dof_encoding = pool;

    return SUANPAN_SUCCESS;
}

void RigidWallPenalty::commit_status() { resistance.reset(); }

void RigidWallPenalty::clear_status() { resistance.reset(); }

void RigidWallPenalty::reset_status() { resistance.reset(); }
