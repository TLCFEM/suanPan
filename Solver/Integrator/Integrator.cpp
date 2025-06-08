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

#include "Integrator.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

Integrator::Integrator(const unsigned T)
    : UniqueTag(T) {}

void Integrator::set_domain(const std::weak_ptr<DomainBase>& D) {
    if(database.lock() != D.lock()) database = D;
}

shared_ptr<DomainBase> Integrator::get_domain() const { return database.lock(); }

int Integrator::initialize() {
    if(nullptr != database.lock()) return SUANPAN_SUCCESS;

    suanpan_error("A valid domain is required.\n");
    return SUANPAN_FAIL;
}

void Integrator::set_time_step_switch(const bool T) { time_step_switch = T; }

/**
 * Some time integration methods (multistep methods) require
 * time step to be constant (for at least some consecutive steps).
 * Call this method in solvers to determine whether it is allowed to change time step.
 */
bool Integrator::allow_to_change_time_step() const { return time_step_switch; }

void Integrator::set_matrix_assembled_switch(const bool T) { matrix_assembled_switch = T; }

bool Integrator::matrix_is_assembled() const { return matrix_assembled_switch; }

bool Integrator::has_corrector() const { return false; }

bool Integrator::time_independent_matrix() const { return true; }

int Integrator::process_load() { return database.lock()->process_load(true); }

/**
 * The main task of this method is to apply constraints (of various forms implemented in various methods).
 * Combinations of different types need to be considered: 1) homogeneous, 2) inhomogeneous, 3) linear, 4) nonlinear.
 * Combinations of different methods need to be considered: 1) penalty, 2) multiplier.
 * On exit, the global stiffness matrix should be updated, the global residual vector should be updated.
 */
int Integrator::process_constraint() {
    const auto D = database.lock();
    auto& W = D->get_factory();

    const auto code = D->process_constraint(true);

    W->set_sushi(W->get_sushi() + W->get_trial_constraint_resistance());

    // some constraints may have stiffness
    D->assemble_constraint_stiffness();

    return code;
}

int Integrator::process_criterion() const { return database.lock()->process_criterion(); }

int Integrator::process_modifier() const { return database.lock()->process_modifier(); }

int Integrator::process_load_resistance() { return database.lock()->process_load(false); }

/**
 * This method is similar to process_constraint(), but it only updates the global residual vector.
 * The global stiffness matrix is not touched as in some solving schemes,
 * the global stiffness matrix is only assembled and factorised once at the beginning.
 * Subsequent iterations do not assemble the global stiffness matrix again and reuse the factorised matrix.
 * In this case, the factorised matrix cannot be modified.
 */
int Integrator::process_constraint_resistance() {
    const auto D = database.lock();
    auto& W = D->get_factory();

    const auto code = D->process_constraint(false);

    W->set_sushi(W->get_sushi() + W->get_trial_constraint_resistance());

    return code;
}

void Integrator::record() const { database.lock()->record(); }

void Integrator::assemble_resistance() {
    const auto D = database.lock();
    auto& W = D->get_factory();
    D->assemble_resistance();
    W->set_sushi(W->get_trial_resistance());
}

/**
 * Assemble the global effective matrix A in AX=B.
 * For FEM applications, it is often a linear combination of stiffness, mass, damping and geometry matrices.
 */
void Integrator::assemble_matrix() {
    const auto D = database.lock();
    auto& W = D->get_factory();
    D->assemble_trial_stiffness();
    D->assemble_trial_geometry();
    if(W->is_nlgeom()) W->get_stiffness() += W->get_geometry();
}

/**
 * Assemble the global residual vector in load-controlled solving schemes.
 */
vec Integrator::get_force_residual() {
    const auto D = database.lock();
    auto& W = D->get_factory();

    vec residual = W->get_trial_load() - W->get_sushi();
    for(const auto I : D->get_restrained_dof()) residual(I) = 0.;
    return residual;
}

/**
 * Assemble the global residual vector in displacement-controlled solving schemes.
 * Apart from the global resistance and external load vectors, the reference load vector shall also be considered.
 */
vec Integrator::get_displacement_residual() {
    const auto D = database.lock();
    auto& W = D->get_factory();

    vec residual = W->get_reference_load() * W->get_trial_load_factor() + W->get_trial_load() - W->get_sushi();
    for(const auto I : D->get_restrained_dof()) residual(I) = 0.;
    return residual;
}

/**
 * Assemble the global residual vector due to nonlinear constraints implemented via the multiplier method.
 */
vec Integrator::get_auxiliary_residual() const {
    auto& W = get_domain()->get_factory();

    return W->get_auxiliary_load() - W->get_auxiliary_resistance();
}

sp_mat Integrator::get_reference_load() { return database.lock()->get_factory()->get_reference_load(); }

const vec& Integrator::get_trial_displacement() const { return database.lock()->get_factory()->get_trial_displacement(); }

void Integrator::update_load() const { database.lock()->update_load(); }

void Integrator::update_constraint() const { database.lock()->update_constraint(); }

void Integrator::update_trial_load_factor(const double lambda) const { update_trial_load_factor(vec{lambda}); }

void Integrator::update_trial_load_factor(const vec& lambda) const {
    auto& W = get_domain()->get_factory();
    W->update_trial_load_factor_by(lambda);
}

void Integrator::update_from_ninja() {
    auto& W = get_domain()->get_factory();
    W->update_trial_displacement_by(W->get_ninja());
}

void Integrator::update_trial_time(const double T) {
    auto& W = get_domain()->get_factory();
    W->update_trial_time(T);
    update_parameter(W->get_incre_time());
}

void Integrator::update_incre_time(const double T) {
    auto& W = get_domain()->get_factory();
    W->update_incre_time(T);
    update_parameter(W->get_incre_time());
}

int Integrator::update_trial_status(const bool detect_trivial) {
    const auto D = get_domain();
    auto& W = D->get_factory();

    return detect_trivial && suanpan::approx_equal(norm(W->get_incre_displacement()), 0.) ? SUANPAN_SUCCESS : D->update_trial_status();
}

int Integrator::correct_trial_status() { return SUANPAN_SUCCESS; }

/**
 * When a new displacement increment is computed, it is added to global displacement vector.
 * At this moment, nodal and elemental quantities are all computed from the previous displacement
 * vector, directly committing the new results causes out-of-sync issue.
 * Some algorithms use predictor-corrector type scheme, which means the converged quantities are
 * different from the committed quantities.
 * This method is in charge of syncing quantities between global and local quantities by updating
 * nodal/elemental quantities using the committed quantities.
 */
int Integrator::sync_status(const bool only_correct) {
    auto handle_force = [&] {
        // process modifiers
        if(SUANPAN_SUCCESS != process_modifier()) return SUANPAN_FAIL;
        // assemble resistance
        assemble_resistance();
        return SUANPAN_SUCCESS;
    };

    // only perform corrector if defined
    if(only_correct) {
        if(!has_corrector()) return SUANPAN_SUCCESS;

        if(SUANPAN_SUCCESS != correct_trial_status()) return SUANPAN_FAIL;

        return handle_force();
    }

    // perform corrector/predictor depending on the algorithm
    if(SUANPAN_SUCCESS != (has_corrector() ? correct_trial_status() : update_trial_status(true))) return SUANPAN_FAIL;

    return handle_force();
}

/**
 * Some algorithms solve a system which differs from the original one.
 * The size of the problem changes thus the computed increment contains additional internal
 * quantities. This method updates internal quantities stored in those integrators.
 */
int Integrator::update_internal(const mat&) { return SUANPAN_SUCCESS; }

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
    if(solve(X, std::move(B)) != SUANPAN_SUCCESS) X.reset();
    return X;
}

mat Integrator::solve(sp_mat&& B) {
    mat X;
    if(solve(X, std::move(B)) != SUANPAN_SUCCESS) X.reset();
    return X;
}

int Integrator::solve(mat& X, const mat& B) { return database.lock()->get_factory()->get_stiffness()->solve(X, B); }

int Integrator::solve(mat& X, const sp_mat& B) { return database.lock()->get_factory()->get_stiffness()->solve(X, B); }

int Integrator::solve(mat& X, mat&& B) { return database.lock()->get_factory()->get_stiffness()->solve(X, std::move(B)); }

int Integrator::solve(mat& X, sp_mat&& B) { return database.lock()->get_factory()->get_stiffness()->solve(X, std::move(B)); }

/**
 * Avoid machine error accumulation.
 * The penalty method can apply homogeneous constraints approximately.
 * The corresponding DoF shall be set to zero after solving the system.
 */
void Integrator::erase_machine_error(vec& ninja) const {
    const auto D = get_domain();
    auto& W = D->get_factory();

    D->erase_machine_error(ninja);
    W->modify_ninja() = ninja.head(W->get_size());
}

void Integrator::stage_and_commit_status() {
    stage_status();
    commit_status();
}

void Integrator::stage_status() const { database.lock()->stage_status(); }

void Integrator::commit_status() { database.lock()->commit_status(); }

void Integrator::clear_status() {
    matrix_assembled_switch = false;
    database.lock()->clear_status();
}

void Integrator::reset_status() { database.lock()->reset_status(); }

/**
 * When time step changes, some parameters may need to be updated.
 */
void Integrator::update_parameter(double) {}

/**
 * When external loads are applied, they can be applied in forms of displacement/velocity/acceleration.
 * The time integration methods, by default, form effective stiffness matrices in displacement domain.
 * That is, in AX=B, A is the effective stiffness matrix and X is the displacement increment.
 * Thus, loads in velocity/acceleration must be converted to displacement.
 * This cannot be done arbitrarily due to compatibility issues.
 * This method takes velocity increment and converts it to TOTAL displacement.
 */
vec Integrator::from_incre_velocity(const vec&, const uvec& encoding) { return zeros(encoding.n_elem); }

/**
 * When external loads are applied, they can be applied in forms of displacement/velocity/acceleration.
 * The time integration methods, by default, form effective stiffness matrices in displacement domain.
 * That is, in AX=B, A is the effective stiffness matrix and X is the displacement increment.
 * Thus, loads in velocity/acceleration must be converted to displacement.
 * This cannot be done arbitrarily due to compatibility issues.
 * This method takes acceleration increment and converts it to TOTAL displacement.
 */
vec Integrator::from_incre_acceleration(const vec&, const uvec& encoding) { return zeros(encoding.n_elem); }

vec Integrator::from_total_velocity(const vec& total_velocity, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    if(AnalysisType::DYNAMICS != W->get_analysis_type()) return zeros(encoding.n_elem);

    return from_incre_velocity(total_velocity - W->get_current_velocity()(encoding), encoding);
}

vec Integrator::from_total_acceleration(const vec& total_acceleration, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

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

bool ImplicitIntegrator::time_independent_matrix() const { return false; }

void ExplicitIntegrator::assemble_resistance() {
    const auto D = get_domain();
    auto& W = D->get_factory();

    auto fa = std::async([&] { D->assemble_resistance(); });
    auto fb = std::async([&] { D->assemble_damping_force(); });
    auto fc = std::async([&] { D->assemble_nonviscous_force(); });
    auto fd = std::async([&] { D->assemble_inertial_force(); });

    fa.get();
    fb.get();
    fc.get();
    fd.get();

    W->set_sushi(W->get_trial_resistance() + W->get_trial_damping_force() + W->get_trial_nonviscous_force() + W->get_trial_inertial_force());
}

void ExplicitIntegrator::assemble_matrix() { get_domain()->assemble_trial_mass(); }

const vec& ExplicitIntegrator::get_trial_displacement() const { return get_domain()->get_factory()->get_trial_acceleration(); }

void ExplicitIntegrator::update_from_ninja() {
    const auto& W = get_domain()->get_factory();
    W->update_trial_acceleration_by(W->get_ninja());
}

int ExplicitIntegrator::solve(mat& X, const mat& B) { return get_domain()->get_factory()->get_mass()->solve(X, B); }

int ExplicitIntegrator::solve(mat& X, const sp_mat& B) { return get_domain()->get_factory()->get_mass()->solve(X, B); }

int ExplicitIntegrator::solve(mat& X, mat&& B) { return get_domain()->get_factory()->get_mass()->solve(X, std::move(B)); }

int ExplicitIntegrator::solve(mat& X, sp_mat&& B) { return get_domain()->get_factory()->get_mass()->solve(X, std::move(B)); }

vec ExplicitIntegrator::from_incre_velocity(const vec&, const uvec&) { throw std::invalid_argument("support velocity cannot be used with explicit integrator"); }

vec ExplicitIntegrator::from_incre_acceleration(const vec& incre_acceleration, const uvec& encoding) { return get_domain()->get_factory()->get_current_acceleration()(encoding) + incre_acceleration; }

vec ExplicitIntegrator::from_total_acceleration(const vec& total_acceleration, const uvec&) { return total_acceleration; }
