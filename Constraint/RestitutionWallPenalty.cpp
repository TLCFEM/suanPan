/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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
#include <Solver/Integrator/Integrator.h>
#include <Step/Step.h>

RestitutionWallPenalty::RestitutionWallPenalty(const unsigned T, const unsigned A, vec&& O, vec&& N, const double RC, const double F, const unsigned NS)
    : RigidWallPenalty(T, A, std::move(O), std::move(N), F, NS)
    , restitution_coefficient(std::max(0., std::min(1., RC))) {}

RestitutionWallPenalty::RestitutionWallPenalty(const unsigned T, const unsigned A, vec&& O, vec&& E1, vec&& E2, const double RC, const double F, const unsigned NS)
    : RigidWallPenalty(T, A, std::move(O), std::move(E1), std::move(E2), F, NS)
    , restitution_coefficient(std::max(0., std::min(1., RC))) {}

int RestitutionWallPenalty::initialize(const shared_ptr<DomainBase>& D) {
    if(AnalysisType::DYNAMICS != D->get_factory()->get_analysis_type()) {
        suanpan_error("Restitution rigid wall constraint can only be applied in dynamic analysis.\n");
        return SUANPAN_FAIL;
    }

    return RigidWallPenalty::initialize(D);
}

int RestitutionWallPenalty::process(const shared_ptr<DomainBase>& D) {
    resistance.reset();
    stiffness.reset();

    suanpan::for_all(D->get_node_pool(), [&](const shared_ptr<Node>& t_node) {
        if(!t_node->validate_dof(ref_dof)) return;
        const vec t_pos = t_node->trial_position(n_dim) - origin;
        if(!edge_a.empty())
            if(const auto projection = dot(t_pos, edge_a); projection > length_a || projection < 0.) return;
        if(!edge_b.empty())
            if(const auto projection = dot(t_pos, edge_b); projection > length_b || projection < 0.) return;
        if(dot(t_pos, outer_norm) > 0.) return;
        node_pool.insert(t_node);
    });

    if(node_pool.empty()) return SUANPAN_SUCCESS;

    auto& W = D->get_factory();
    auto& G = D->get_current_step()->get_integrator();

    const auto factor = alpha * pow(W->get_incre_time(), -2.);

    std::vector<uword> pool;
    pool.reserve(n_dim * node_pool.size());

    resistance.zeros(W->get_size());
    auto counter = 0llu;
    for(const auto& I : node_pool) {
        const auto c_vel = I->get_current_velocity(n_dim);
        if(dot(c_vel, outer_norm) > 0.) continue;
        auto& t_dof = I->get_reordered_dof();
        const auto t_vel = I->get_trial_velocity(n_dim);
        const vec diff_disp = I->get_trial_displacement(n_dim) - G->from_total_velocity(t_vel - dot(t_vel + restitution_coefficient * c_vel, outer_norm) * outer_norm, t_dof.head(n_dim));
        const auto next_counter = counter + n_dim;
        stiffness.resize(next_counter, next_counter);
        stiffness.submat(counter, counter, size(n_dim, n_dim)) = factor * outer_norm * outer_norm.t();
        for(auto J = 0llu; J < n_dim; ++J) {
            pool.emplace_back(t_dof(J));
            resistance(t_dof(J)) += factor * diff_disp(J);
        }
        counter = next_counter;
    }

    dof_encoding = pool;

    return SUANPAN_SUCCESS;
}

void RestitutionWallPenalty::stage(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    auto trial_acceleration = W->modify_trial_acceleration();
    for(const auto& I : node_pool) {
        auto t_acceleration = I->get_trial_acceleration();
        t_acceleration.head(n_dim) = t_acceleration.head(n_dim) - dot(I->get_incre_acceleration(n_dim), outer_norm) * outer_norm;
        trial_acceleration(I->get_reordered_dof()) = I->update_trial_acceleration(std::move(t_acceleration));
    }

    W->update_trial_acceleration(trial_acceleration);
}

void RestitutionWallPenalty::commit_status() { node_pool.clear(); }

void RestitutionWallPenalty::clear_status() { node_pool.clear(); }

void RestitutionWallPenalty::reset_status() { node_pool.clear(); }

RestitutionWallPenalty1D::RestitutionWallPenalty1D(const unsigned T, const unsigned A, const double O, const double N, const double RC, const double F)
    : RestitutionWallPenalty(T, A, {O}, {N}, RC, F, 1) {}

RestitutionWallPenalty2D::RestitutionWallPenalty2D(const unsigned T, const unsigned A, vec2&& O, vec2&& N, const double RC, const double F)
    : RestitutionWallPenalty(T, A, std::move(O), std::move(N), RC, F, 2) {}

RestitutionWallPenalty2D::RestitutionWallPenalty2D(const unsigned T, const unsigned A, vec2&& O, vec3&& E1, const double RC, const double F)
    : RestitutionWallPenalty(T, A, std::move(O), std::move(E1), vec{0., 0., 1.}, RC, F, 2) {
    access::rw(outer_norm).resize(2);
    access::rw(edge_a).resize(2);
    access::rw(edge_b).reset();
}

RestitutionWallPenalty3D::RestitutionWallPenalty3D(const unsigned T, const unsigned A, vec3&& O, vec3&& N, const double RC, const double F)
    : RestitutionWallPenalty(T, A, std::move(O), std::move(N), RC, F, 3) {}

RestitutionWallPenalty3D::RestitutionWallPenalty3D(const unsigned T, const unsigned A, vec3&& O, vec3&& E1, vec3&& E2, const double RC, const double F)
    : RestitutionWallPenalty(T, A, std::move(O), std::move(E1), std::move(E2), RC, F, 3) {}
