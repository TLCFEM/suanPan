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

#include "Integrator.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

Integrator::Integrator(const unsigned T)
    : Tag(T) { suanpan_debug("Integrator %u ctor() called.\n", T); }

Integrator::~Integrator() { suanpan_debug("Integrator %u dtor() called.\n", get_tag()); }

void Integrator::set_domain(const weak_ptr<DomainBase>& D) { if(database.lock() != D.lock()) database = D; }

const weak_ptr<DomainBase>& Integrator::get_domain() const { return database; }

int Integrator::initialize() {
    if(nullptr != database.lock()) return SUANPAN_SUCCESS;

    suanpan_error("initialize() needs a valid domain.\n");
    return SUANPAN_FAIL;
}

void Integrator::set_time_step_switch(const bool T) { time_step_switch = T; }

/**
 * Some time integration methods (multistep methods) require
 * time step to be constant (for at least some consecutive steps).
 * Call this method in solvers to determine whether it is allowed to change time step.
 */
bool Integrator::allow_to_change_time_step() const { return time_step_switch; }

int Integrator::process_load() { return database.lock()->process_load(true); }

/**
 * The main task of this method is to apply constraints (of various forms implemented in various methods).
 * Combinations of different types need to be considered: 1) homogeneous, 2) inhomogeneous, 3) linear, 4) nonlinear.
 * Combinations of different methods need to be considered: 1) penalty, 2) multiplier.
 * On exit, the global stiffness matrix should be updated, the global residual vector should be updated.
 */
int Integrator::process_constraint() {
    const auto& D = database.lock();
    const auto& W = D->get_factory();

    const auto code = D->process_constraint(true);

    W->set_sushi(W->get_sushi() + W->get_trial_constraint_resistance());

    // some constraints may have stiffness
    D->assemble_constraint_stiffness();

    return code;
}

int Integrator::process_criterion() { return database.lock()->process_criterion(); }

int Integrator::process_modifier() { return database.lock()->process_modifier(); }

int Integrator::process_load_resistance() { return SUANPAN_SUCCESS; }

/**
 * This method is similar to process_constraint(), but it only updates the global residual vector.
 * The global stiffness matrix is not touched as in some solving schemes,
 * the global stiffness matrix is only assembled and factorised once at the beginning.
 * Subsequent iterations do not assemble the global stiffness matrix again and reuse the factorised matrix.
 * In this case, the factorised matrix cannot be modified.
 */
int Integrator::process_constraint_resistance() {
    const auto& D = database.lock();
    const auto& W = D->get_factory();

    const auto code = D->process_constraint(false);

    W->set_sushi(W->get_sushi() + W->get_trial_constraint_resistance());

    return code;
}

void Integrator::record() const { database.lock()->record(); }

void Integrator::assemble_resistance() {
    const auto& D = database.lock();
    const auto& W = D->get_factory();
    D->assemble_resistance();
    W->set_sushi(W->get_trial_resistance());
}

/**
 * Assemble the global effective matrix A in AX=B.
 * For FEM applications, it is often a linear combination of stiffness, mass, damping and geometry matrices.
 */
void Integrator::assemble_matrix() {
    const auto& D = database.lock();
    auto& W = D->get_factory();
    D->assemble_trial_stiffness();
    D->assemble_trial_geometry();
    W->get_stiffness() += W->get_geometry();
}

/**
 * Assemble the global residual vector in load-controlled solving schemes.
 */
vec Integrator::get_force_residual() {
    const auto& D = database.lock();
    const auto& W = D->get_factory();

    vec residual = W->get_trial_load() - W->get_sushi();

    for(auto& I : D->get_restrained_dof()) residual(I) = 0.;

    return residual;
}

/**
 * Assemble the global residual vector in displacement-controlled solving schemes.
 * Apart from the global resistance and external load vectors, the reference load vector shall also be considered.
 */
vec Integrator::get_displacement_residual() {
    const auto& D = database.lock();
    const auto& W = D->get_factory();

    vec residual = W->get_reference_load() * W->get_trial_load_factor() + W->get_trial_load() - W->get_sushi();

    for(auto& I : D->get_restrained_dof()) residual(I) = 0.;

    return residual;
}

/**
 * Assemble the global residual vector due to nonlinear constraints implemented via the multiplier method.
 */
vec Integrator::get_auxiliary_residual() {
    const auto& W = get_domain().lock()->get_factory();

    return W->get_auxiliary_load() - W->get_auxiliary_resistance();
}

sp_mat Integrator::get_reference_load() { return database.lock()->get_factory()->get_reference_load(); }

sp_mat Integrator::get_auxiliary_stiffness() { return database.lock()->get_factory()->get_auxiliary_stiffness(); }

void Integrator::update_load() { database.lock()->update_load(); }

void Integrator::update_constraint() { database.lock()->update_constraint(); }

void Integrator::update_trial_load_factor(const double lambda) { update_trial_load_factor(vec{lambda}); }

void Integrator::update_trial_load_factor(const vec& lambda) {
    const auto& W = get_domain().lock()->get_factory();
    W->update_trial_load_factor(W->get_trial_load_factor() + lambda);
}

void Integrator::update_trial_displacement(const vec& ninja) {
    const auto& W = get_domain().lock()->get_factory();
    W->update_trial_displacement(W->get_trial_displacement() + ninja);
}

void Integrator::update_trial_time(const double T) {
    const auto& W = get_domain().lock()->get_factory();
    W->update_trial_time(T);
    update_parameter(W->get_incre_time());
}

void Integrator::update_incre_time(const double T) {
    const auto& W = get_domain().lock()->get_factory();
    W->update_incre_time(T);
    update_parameter(W->get_incre_time());
}

int Integrator::update_trial_status() { return database.lock()->update_trial_status(); }

int Integrator::update_incre_status() { return database.lock()->update_incre_status(); }

/**
 * Must change ninja to the real displacement increment.
 */
int Integrator::update_internal(const mat&) { return 0; }

mat Integrator::solve(const mat& B) {
    mat X;
    if(solve(X, B) != SUANPAN_SUCCESS) X.reset();
    return X;
}

mat Integrator::solve(const sp_mat& B) {
    mat X;
    if(solve(X, B) != SUANPAN_SUCCESS) X.reset();
    return X;
}

mat Integrator::solve(mat&& B) {
    mat X;
    if(solve(X, std::forward<mat>(B)) != SUANPAN_SUCCESS) X.reset();
    return X;
}

mat Integrator::solve(sp_mat&& B) {
    mat X;
    if(solve(X, std::forward<sp_mat>(B)) != SUANPAN_SUCCESS) X.reset();
    return X;
}

int Integrator::solve(mat& X, const mat& B) { return database.lock()->get_factory()->get_stiffness()->solve(X, B); }

int Integrator::solve(mat& X, const sp_mat& B) { return database.lock()->get_factory()->get_stiffness()->solve(X, B); }

int Integrator::solve(mat& X, mat&& B) { return database.lock()->get_factory()->get_stiffness()->solve(X, std::forward<mat>(B)); }

int Integrator::solve(mat& X, sp_mat&& B) { return database.lock()->get_factory()->get_stiffness()->solve(X, std::forward<sp_mat>(B)); }

/**
 * Avoid machine error accumulation.
 * The penalty method can apply homogeneous constraints approximately.
 * The corresponding DoF shall be set to zero after solving the system.
 */
void Integrator::erase_machine_error() const { database.lock()->erase_machine_error(); }

void Integrator::stage_and_commit_status() {
    stage_status();
    commit_status();
}

void Integrator::stage_status() { database.lock()->stage_status(); }

void Integrator::commit_status() {
    database.lock()->commit_status();
    update_compatibility();
}

void Integrator::clear_status() {
    database.lock()->clear_status();
    update_compatibility();
}

void Integrator::reset_status() {
    database.lock()->reset_status();
    update_compatibility();
}

/**
 * When tim step changes, some parameters may need to be updated.
 */
void Integrator::update_parameter(double) {}

/**
 * Make sure that the trial displacement/velocity/acceleration are consistent with each other.
 * When starting a new trial state, the trial displacement is identical to the current displacement.
 * This essentially means that the displacement increment is zero.
 * To have such a trial state with the given time step, the trial velocity and acceleration shall be
 * updated to be compatible with the trial displacement.
 *
 * On exit, trial velocity and acceleration should be computed from current/trial displacement.
 */
void Integrator::update_compatibility() const {}

/**
 * When external loads are applied, they can be applied in forms of displacement/velocity/acceleration.
 * The time integration methods, by default, form effective stiffness matrices in displacement domain.
 * That is, in AX=B, A is the effective stiffness matrix and X is the displacement increment.
 * Thus, loads in velocity/acceleration must be converted to displacement.
 * This cannot be done arbitrarily due to compatibility issues.
 * This method takes velocity increment and converts it to displacement increment.
 */
vec Integrator::from_incre_velocity(const vec&, const uvec& encoding) { return zeros(encoding.n_elem); }

/**
 * When external loads are applied, they can be applied in forms of displacement/velocity/acceleration.
 * The time integration methods, by default, form effective stiffness matrices in displacement domain.
 * That is, in AX=B, A is the effective stiffness matrix and X is the displacement increment.
 * Thus, loads in velocity/acceleration must be converted to displacement.
 * This cannot be done arbitrarily due to compatibility issues.
 * This method takes acceleration increment and converts it to displacement increment.
 */
vec Integrator::from_incre_acceleration(const vec&, const uvec& encoding) { return zeros(encoding.n_elem); }

vec Integrator::from_total_velocity(const vec& total_velocity, const uvec& encoding) {
    const auto& W = get_domain().lock()->get_factory();

    if(AnalysisType::DYNAMICS != W->get_analysis_type()) return zeros(encoding.n_elem);

    return from_incre_velocity(total_velocity - W->get_current_velocity()(encoding), encoding);
}

vec Integrator::from_total_acceleration(const vec& total_acceleration, const uvec& encoding) {
    const auto& W = get_domain().lock()->get_factory();

    if(AnalysisType::DYNAMICS != W->get_analysis_type()) return zeros(encoding.n_elem);

    return from_incre_acceleration(total_acceleration - W->get_current_acceleration()(encoding), encoding);
}

/**
 * A simplified version similar to `from_incre_velocity(const vec&, const uvec&)`.
 * It assumes all DoFs share the same magnitude.
 */
vec Integrator::from_incre_velocity(const double magnitude, const uvec& encoding) { return from_incre_velocity(vec(encoding.n_elem, fill::value(magnitude)), encoding); }

/**
 * A simplified version similar to `from_incre_acceleration(const vec&, const uvec&)`.
 * It assumes all DoFs share the same magnitude.
 */
vec Integrator::from_incre_acceleration(const double magnitude, const uvec& encoding) { return from_incre_acceleration(vec(encoding.n_elem, fill::value(magnitude)), encoding); }

vec Integrator::from_total_velocity(const double magnitude, const uvec& encoding) { return from_total_velocity(vec(encoding.n_elem, fill::value(magnitude)), encoding); }

vec Integrator::from_total_acceleration(const double magnitude, const uvec& encoding) { return from_total_acceleration(vec(encoding.n_elem, fill::value(magnitude)), encoding); }
