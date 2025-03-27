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

#include "LeeNewmarkBase.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

int LeeNewmarkBase::erase_top_left_block() const {
    auto& t_triplet = stiffness->triplet_mat;

    uword *ptr_a, *ptr_b;

    if(t_triplet.is_csc_sorted()) {
        ptr_a = t_triplet.col_mem();
        ptr_b = t_triplet.row_mem();
    }
    else if(t_triplet.is_csr_sorted()) {
        ptr_a = t_triplet.row_mem();
        ptr_b = t_triplet.col_mem();
    }
    else {
        suanpan_error("The system is not sorted, please file a bug report.\n");
        return SUANPAN_FAIL;
    }

    const auto& val = t_triplet.val_mem();

    for(uword I = 0; I < t_triplet.n_elem; ++I) {
        // quit if current column/row is beyond the original size of matrix
        if(ptr_a[I] >= n_block) break;
        // erase existing entries if fall in intact stiffness matrix
        if(ptr_b[I] < n_block) val[I] = 0.;
    }

    return SUANPAN_SUCCESS;
}

LeeNewmarkBase::LeeNewmarkBase(const unsigned T, const double A, const double B, const StiffnessType ST)
    : Newmark(T, A, B)
    , n_block(0)
    , stiffness_type(ST) {}

int LeeNewmarkBase::initialize() {
    if(Newmark::initialize() != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    factory = get_domain()->get_factory();

    access::rw(n_block) = factory->get_size();

    const auto n_size = get_total_size();

    trial_internal = current_internal = residual.zeros(n_size);

    if(factory->contain_solver_type(SolverType::LIS)) stiffness = make_unique<SparseMatLis<double>>(n_size, n_size);
#ifndef SUANPAN_DISTRIBUTED
    else if(factory->contain_solver_type(SolverType::MUMPS)) stiffness = make_unique<SparseMatMUMPS<double>>(n_size, n_size);
#endif
#ifdef SUANPAN_MKL
    else if(factory->contain_solver_type(SolverType::PARDISO)) stiffness = make_unique<SparseMatPARDISO<double>>(n_size, n_size);
    else if(factory->contain_solver_type(SolverType::FGMRES)) stiffness = make_unique<SparseMatFGMRES<double>>(n_size, n_size);
#endif
#ifdef SUANPAN_CUDA
    else if(factory->contain_solver_type(SolverType::CUDA)) stiffness = make_unique<SparseMatCUDA<double>>(n_size, n_size);
#endif
    else stiffness = make_unique<SparseMatSuperLU<double>>(n_size, n_size);

    return SUANPAN_SUCCESS;
}

int LeeNewmarkBase::update_internal(const mat& t_internal) {
    trial_internal += t_internal;

    return SUANPAN_SUCCESS;
}

int LeeNewmarkBase::solve(mat& X, const mat& B) {
    stiffness->set_solver_setting(factory->get_solver_setting());
    return stiffness->solve(X, resize(B, stiffness->n_rows, B.n_cols));
}

int LeeNewmarkBase::solve(mat& X, const sp_mat& B) {
    stiffness->set_solver_setting(factory->get_solver_setting());
    return stiffness->solve(X, resize(B, stiffness->n_rows, B.n_cols));
}

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
