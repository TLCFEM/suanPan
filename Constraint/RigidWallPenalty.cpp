/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#include <Domain/Factory.hpp>

RigidWallPenalty::RigidWallPenalty(const unsigned T, const unsigned A, vec&& O, vec&& N, const double F, const unsigned NS)
    : Constraint(T, A, setup(NS), {}, 0)
    , n_dim(NS)
    , alpha(F)
    , origin(std::move(O))
    , outer_norm(normalise(N)) {}

RigidWallPenalty::RigidWallPenalty(const unsigned T, const unsigned A, vec&& O, vec&& E1, vec&& E2, const double F, const unsigned NS)
    : Constraint(T, A, setup(NS), {}, 0)
    , n_dim(NS)
    , alpha(F)             // penalty factor
    , edge_a(E1)           // 3D vector
    , edge_b(E2)           // 3D vector
    , origin(std::move(O)) // 1D, 2D and 3D vectors
    , outer_norm(normalise(cross(edge_a, edge_b)))
    , length_a(norm(edge_a))
    , length_b(norm(edge_b)) {}

int RigidWallPenalty::initialize(const shared_ptr<DomainBase>& D) {
    if(SUANPAN_SUCCESS != Constraint::initialize(D)) return SUANPAN_FAIL;

    return validate_node(D) ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

int RigidWallPenalty::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    const auto factor = alpha * pow(W->get_incre_time(), -2.);

    stiffness.reset();
    resistance.zeros(W->get_size());
    std::vector<uword> pool;

    uword counter{0};
    for(const auto& I : D->get_node_pool()) {
        const vec t_pos = I->trial_position(n_dim) - origin;
        if(!edge_a.empty())
            if(const auto projection = dot(t_pos, edge_a); projection > length_a || projection < 0.) continue;
        if(!edge_b.empty())
            if(const auto projection = dot(t_pos, edge_b); projection > length_b || projection < 0.) continue;
        const auto t_pen = dot(t_pos, outer_norm);
        if(t_pen > datum::eps) continue;
        const auto next_counter = counter + n_dim;
        stiffness.resize(next_counter, next_counter);
        stiffness.submat(counter, counter, size(n_dim, n_dim)) = factor * outer_norm * outer_norm.t();
        auto& t_dof = I->get_reordered_dof();
        for(uword J{0}; J < n_dim; ++J) {
            pool.emplace_back(t_dof(J));
            resistance(t_dof(J)) += factor * t_pen * outer_norm(J);
        }
        counter = next_counter;
    }

    target_dof = pool;

    return SUANPAN_SUCCESS;
}

void RigidWallPenalty::commit_status() { resistance.reset(); }

void RigidWallPenalty::clear_status() { resistance.reset(); }

void RigidWallPenalty::reset_status() { resistance.reset(); }

RigidWallPenalty1D::RigidWallPenalty1D(const unsigned T, const unsigned A, const double O, const double N, const double F)
    : RigidWallPenalty(T, A, {O}, {N}, F, 1) {}

RigidWallPenalty2D::RigidWallPenalty2D(const unsigned T, const unsigned A, vec2&& O, vec2&& N, const double F)
    : RigidWallPenalty(T, A, std::move(O), std::move(N), F, 2) {}

RigidWallPenalty2D::RigidWallPenalty2D(const unsigned T, const unsigned A, vec2&& O, vec3&& E1, const double F)
    : RigidWallPenalty(T, A, std::move(O), std::move(E1), vec{0., 0., 1.}, F, 2) {
    access::rw(outer_norm).resize(2);
    access::rw(edge_a).resize(2);
    access::rw(edge_b).reset();
}

RigidWallPenalty3D::RigidWallPenalty3D(const unsigned T, const unsigned A, vec3&& O, vec3&& N, const double F)
    : RigidWallPenalty(T, A, std::move(O), std::move(N), F, 3) {}

RigidWallPenalty3D::RigidWallPenalty3D(const unsigned T, const unsigned A, vec3&& O, vec3&& E1, vec3&& E2, const double F)
    : RigidWallPenalty(T, A, std::move(O), std::move(E1), std::move(E2), F, 3) {}
