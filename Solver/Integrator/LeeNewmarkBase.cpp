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

#include "LeeNewmarkBase.h"
#include <Domain/DomainBase.h>
#include <Domain/FactoryHelper.hpp>

LeeNewmarkBase::LeeNewmarkBase(const unsigned T, const double A, const double B, const StiffnessType ST)
    : Newmark(T, A, B)
    , n_block(0)
    , stiffness_type(ST) {}

int LeeNewmarkBase::initialize() {
    if(Newmark::initialize() != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    factory = get_domain().lock()->get_factory();

    access::rw(n_block) = factory->get_size();

    const auto n_size = get_total_size();

    trial_internal = current_internal = residual.zeros(n_size);

    if(SolverType::MUMPS == factory->get_solver_type()) stiffness = make_unique<SparseMatMUMPS<double>>(n_size, n_size);
#ifdef SUANPAN_MKL
    else if(SolverType::PARDISO == factory->get_solver_type()) stiffness = make_unique<SparseMatPARDISO<double>>(n_size, n_size);
    else if(SolverType::FGMRES == factory->get_solver_type()) stiffness = make_unique<SparseMatFGMRES<double>>(n_size, n_size);
#endif
#ifdef SUANPAN_CUDA
    else if(SolverType::CUDA == factory->get_solver_type()) stiffness = make_unique<SparseMatCUDA<double>>(n_size, n_size);
#endif
    else stiffness = make_unique<SparseMatSuperLU<double>>(n_size, n_size);

    stiffness->set_solver_setting(factory->get_solver_setting());

    return SUANPAN_SUCCESS;
}

int LeeNewmarkBase::update_internal(const mat& t_internal) {
    trial_internal += t_internal;

    // order matters
    // cannot resize before assignment
    get_ninja(factory).resize(n_block);

    return SUANPAN_SUCCESS;
}

int LeeNewmarkBase::solve(mat& X, const mat& B) { return stiffness->solve(X, resize(B, stiffness->n_rows, B.n_cols)); }

int LeeNewmarkBase::solve(mat& X, const sp_mat& B) { return stiffness->solve(X, resize(B, stiffness->n_rows, B.n_cols)); }

int LeeNewmarkBase::solve(mat& X, mat&& B) { return solve(X, B); }

int LeeNewmarkBase::solve(mat& X, sp_mat&& B) { return solve(X, B); }

vec LeeNewmarkBase::get_force_residual() {
    residual.head_rows(n_block) = Newmark::get_force_residual();

    return residual;
}

vec LeeNewmarkBase::get_displacement_residual() {
    residual.head_rows(n_block) = Newmark::get_displacement_residual();

    return residual;
}

void LeeNewmarkBase::commit_status() {
    current_internal = trial_internal;

    first_iteration = true;

    Newmark::commit_status();
}

void LeeNewmarkBase::clear_status() {
    current_internal = trial_internal.zeros();

    first_iteration = true;

    Newmark::clear_status();
}

void LeeNewmarkBase::reset_status() {
    trial_internal = current_internal;

    first_iteration = true;

    Newmark::reset_status();
}
