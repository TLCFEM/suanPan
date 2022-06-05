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

#include "RestitutionWallPenalty.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Solver/Integrator/Integrator.h>
#include <Step/Step.h>

RestitutionWallPenalty::RestitutionWallPenalty(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& N, const double RC, const double F)
    : RigidWallPenalty(T, S, A, std::forward<vec>(O), std::forward<vec>(N), F)
    , restitution_coefficient(std::max(0., std::min(1., RC))) {}

RestitutionWallPenalty::RestitutionWallPenalty(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& E1, vec&& E2, const double RC, const double F)
    : RigidWallPenalty(T, S, A, std::forward<vec>(O), std::forward<vec>(E1), std::forward<vec>(E2), F)
    , restitution_coefficient(std::max(0., std::min(1., RC))) {}

int RestitutionWallPenalty::process(const shared_ptr<DomainBase>& D) {
    resistance.reset();
    stiffness.reset();

    const auto& t_node_pool = D->get_node_pool();
    suanpan_for_each(t_node_pool.cbegin(), t_node_pool.cend(), [&](const shared_ptr<Node>& t_node) {
        auto& t_coor = t_node->get_coordinate();
        auto& t_disp = t_node->get_trial_displacement();
        const auto t_size = std::min(outer_norm.n_elem, std::min(t_coor.n_elem, t_disp.n_elem));
        vec t_pos = -origin;
        for(auto J = 0llu; J < t_size; ++J) t_pos(J) += t_coor(J) + t_disp(J);
        if(!edge_a.empty() && dot(t_pos, edge_a) > length_a || !edge_b.empty() && dot(t_pos, edge_b) > length_b || dot(t_pos, outer_norm) > datum::eps) return;
        node_pool.insert(t_node);
        D->register_node_to_reset_acceleration(t_node->get_tag());
    });

    if(node_pool.empty()) return SUANPAN_SUCCESS;

    auto& W = D->get_factory();
    auto& G = D->get_current_step()->get_integrator();

    const auto factor = alpha * pow(W->get_incre_time(), -2.);

    vector<uword> pool;
    pool.reserve(3llu * node_pool.size());

    resistance.zeros(W->get_size());
    for(const auto& I : node_pool) {
        auto& t_dof = I->get_reordered_dof();
        auto& t_vel = I->get_current_velocity();
        const auto t_size = t_vel.n_elem;
        const vec diff_vel = (-1. - restitution_coefficient) * dot(t_vel, outer_norm.head(t_size)) * outer_norm.head(t_size);
        const vec diff_disp = I->get_trial_displacement() - G->from_incre_velocity(diff_vel, t_dof);
        for(auto J = 0llu; J < t_size; ++J) {
            pool.emplace_back(t_dof(J));
            resistance(t_dof(J)) += factor * diff_disp(J);
        }
    }

    dof_encoding = pool;

    stiffness = speye(dof_encoding.n_elem, dof_encoding.n_elem) * factor;

    return SUANPAN_SUCCESS;
}

void RestitutionWallPenalty::commit_status() { node_pool.clear(); }

void RestitutionWallPenalty::clear_status() { node_pool.clear(); }

void RestitutionWallPenalty::reset_status() { node_pool.clear(); }
