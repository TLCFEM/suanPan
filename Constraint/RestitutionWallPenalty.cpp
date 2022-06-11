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
#include <Domain/FactoryHelper.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Step/Step.h>

RestitutionWallPenalty::RestitutionWallPenalty(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& N, const double RC, const double F, const unsigned NS)
    : RigidWallPenalty(T, S, A, std::forward<vec>(O), std::forward<vec>(N), F, NS)
    , restitution_coefficient(std::max(0., std::min(1., RC))) {}

RestitutionWallPenalty::RestitutionWallPenalty(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& E1, vec&& E2, const double RC, const double F, const unsigned NS)
    : RigidWallPenalty(T, S, A, std::forward<vec>(O), std::forward<vec>(E1), std::forward<vec>(E2), F, NS)
    , restitution_coefficient(std::max(0., std::min(1., RC))) {}

int RestitutionWallPenalty::initialize(const shared_ptr<DomainBase>& D) {
    if(AnalysisType::DYNAMICS != D->get_factory()->get_analysis_type()) {
        suanpan_error("restitution rigid wall constraint can only be applied in dynamic analysis.\n");
        return SUANPAN_FAIL;
    }

    return RigidWallPenalty::initialize(D);
}

int RestitutionWallPenalty::process(const shared_ptr<DomainBase>& D) {
    resistance.reset();
    stiffness.reset();

    suanpan::for_all(D->get_node_pool(), [&](const shared_ptr<Node>& t_node) {
        if(!checker_handler(t_node)) return;
        const vec t_pos = trial_position_handler(t_node) - origin;
        if(!edge_a.empty()) if(const auto projection = dot(t_pos, edge_a); projection > length_a || projection < 0.) return;
        if(!edge_b.empty()) if(const auto projection = dot(t_pos, edge_b); projection > length_b || projection < 0.) return;
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
        if(dot(current_velocity_handler(I), outer_norm) > 0.) continue;
        const auto c_vel = current_velocity_handler(I);
        if(dot(c_vel, outer_norm) > 0.) continue;
        auto& t_dof = I->get_reordered_dof();
        const auto t_vel = trial_velocity_handler(I);
        const vec diff_disp = trial_displacement_handler(I) - G->from_total_velocity(t_vel - dot(t_vel + restitution_coefficient * c_vel, outer_norm) * outer_norm, t_dof.head(n_dim));
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

int RestitutionWallPenalty::stage(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    auto trial_acceleration = get_trial_acceleration(W);
    for(const auto& I : node_pool) {
        auto t_acceleration = I->get_trial_acceleration();
        t_acceleration.head(n_dim) = trial_acceleration_handler(I) - dot(incre_acceleration_handler(I), outer_norm) * outer_norm;
        I->update_trial_acceleration(t_acceleration);
        trial_acceleration(I->get_reordered_dof()) = t_acceleration;
    }

    W->update_trial_acceleration(trial_acceleration);

    return SUANPAN_SUCCESS;
}

void RestitutionWallPenalty::commit_status() { node_pool.clear(); }

void RestitutionWallPenalty::clear_status() { node_pool.clear(); }

void RestitutionWallPenalty::reset_status() { node_pool.clear(); }

RestitutionWallPenalty1D::RestitutionWallPenalty1D(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& N, const double RC, const double F)
    : RestitutionWallPenalty(T, S, A, resize(O, 1, 1), resize(N, 1, 1), RC, F, 1) { set_handler<DOF::U1>(); }

RestitutionWallPenalty2D::RestitutionWallPenalty2D(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& N, const double RC, const double F)
    : RestitutionWallPenalty(T, S, A, resize(O, 2, 1), resize(N, 2, 1), RC, F, 2) { set_handler<DOF::U1, DOF::U2>(); }

RestitutionWallPenalty2D::RestitutionWallPenalty2D(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& E1, vec&& E2, const double RC, const double F)
    : RestitutionWallPenalty(T, S, A, resize(O, 2, 1), resize(E1, 3, 1), resize(E2, 3, 1), RC, F, 2) {
    set_handler<DOF::U1, DOF::U2>();
    access::rw(outer_norm).resize(2);
    access::rw(edge_a).resize(2);
    access::rw(edge_b).reset();
}

RestitutionWallPenalty3D::RestitutionWallPenalty3D(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& N, const double RC, const double F)
    : RestitutionWallPenalty(T, S, A, resize(O, 3, 1), resize(N, 3, 1), RC, F, 3) { set_handler<DOF::U1, DOF::U2, DOF::U3>(); }

RestitutionWallPenalty3D::RestitutionWallPenalty3D(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& E1, vec&& E2, const double RC, const double F)
    : RestitutionWallPenalty(T, S, A, resize(O, 3, 1), resize(E1, 3, 1), resize(E2, 3, 1), RC, F, 3) { set_handler<DOF::U1, DOF::U2, DOF::U3>(); }
