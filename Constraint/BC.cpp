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

#include "BC.h"

#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Step/Step.h>

double PenaltyBC::multiplier = 1E8;

PenaltyBC::PenaltyBC(const unsigned T, uvec&& N, std::vector<Node::DOF>&& D)
    : Constraint(T, 0, {}, std::move(D), 0) { target_node = std::move(N); }

int PenaltyBC::initialize(const shared_ptr<DomainBase>& D) {
    if(SUANPAN_SUCCESS != Constraint::initialize(D)) return SUANPAN_FAIL;

    target_dof = collect_node_dof(D);

    return SUANPAN_SUCCESS;
}

/**
 * \brief Apply the BC to the system using penalty method.
 * It effectively adds a diagonal matrix to the global stiffness matrix.
 */
int PenaltyBC::process(const shared_ptr<DomainBase>& D) {
    if(Integrator::Type::Explicit == D->get_current_step()->get_integrator()->type()) return process_resistance(D);

    auto& W = D->get_factory();

    double reference{0.};
    if(AnalysisType::EIGEN == W->get_analysis_type()) reference = W->get_stiffness()->max();
    else {
        const auto incre_time = W->get_incre_time();
        if(auto& t_mat = W->get_stiffness()) reference += t_mat->max();
        if(auto& t_mat = W->get_geometry()) reference += t_mat->max();
        if(auto& t_mat = W->get_damping()) reference += t_mat->max() / incre_time;
        if(auto& t_mat = W->get_nonviscous()) reference += t_mat->max() / incre_time;
        if(auto& t_mat = W->get_mass()) reference += t_mat->max() / incre_time / incre_time;
    }

    stiffness.zeros(target_dof.n_elem, target_dof.n_elem).diag().fill(multiplier * reference);

    return process_resistance(D);
}

int PenaltyBC::process_resistance(const shared_ptr<DomainBase>& D) {
    // this ensures all restrained DoFs have zero displacement in modified Newton methods
    D->insert_constrained_dof(target_dof);

    return SUANPAN_SUCCESS;
}

/**
 * \brief Apply the BC to the system using Lagrangian multiplier method.
 * It directly modifies the global stiffness matrix thus requires a mutex.
 * Other global matrices shall also be modified to ensure the solution is trivial on target DoFs.
 */
int MultiplierBC::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    if(Integrator::Type::Explicit == D->get_current_step()->get_integrator()->type()) {
        if(auto& t_mass = W->get_mass()) {
            std::scoped_lock lock(W->get_mass_mutex());
            t_mass->unify(target_dof);
        }
    }
    else {
        if(auto& t_stiff = W->get_stiffness()) {
            std::scoped_lock lock(W->get_stiffness_mutex());
            t_stiff->unify(target_dof);
        }
        if(auto& t_mass = W->get_mass()) {
            std::scoped_lock lock(W->get_mass_mutex());
            t_mass->nullify(target_dof);
        }
        if(auto& t_damping = W->get_damping()) {
            std::scoped_lock lock(W->get_damping_mutex());
            t_damping->nullify(target_dof);
        }
        if(auto& t_nonviscous = W->get_nonviscous()) {
            std::scoped_lock lock(W->get_nonviscous_mutex());
            t_nonviscous->nullify(target_dof);
        }
        if(auto& t_geometry = W->get_geometry()) {
            std::scoped_lock lock(W->get_geometry_mutex());
            t_geometry->nullify(target_dof);
        }
    }

    return process_resistance(D);
}

GroupPenaltyBC::GroupPenaltyBC(const unsigned T, uvec&& N, std::vector<Node::DOF>&& D)
    : GroupModifier(std::move(N))
    , MultiplierBC(T, {}, std::move(D)) {}

int GroupPenaltyBC::initialize(const shared_ptr<DomainBase>& D) {
    target_node = update_object_tag(D);

    return MultiplierBC::initialize(D);
}

int GroupPenaltyBC::process(const shared_ptr<DomainBase>& D) {
    return PenaltyBC::process(D); // NOLINT(bugprone-parent-virtual-call)
}

int GroupMultiplierBC::process(const shared_ptr<DomainBase>& D) {
    return MultiplierBC::process(D); // NOLINT(bugprone-parent-virtual-call)
}

void set_constraint_multiplier(const double M) { PenaltyBC::multiplier = M; }
