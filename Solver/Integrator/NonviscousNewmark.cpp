/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include "NonviscousNewmark.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

mat NonviscousNewmark::get_residual() const {
    auto& W = get_domain()->get_factory();

    return trial_damping - current_damping * diagmat(C5 / (C5 + s)) - W->get_trial_velocity() * (m / (C5 + s)).t();
}

NonviscousNewmark::NonviscousNewmark(const unsigned T, const double A, const double B, vec&& M, vec&& S)
    : Newmark(T, A, B)
    , m(std::forward<vec>(M))
    , s(std::forward<vec>(S)) {}

int NonviscousNewmark::initialize() {
    auto& W = get_domain()->get_factory();

    trial_damping = current_damping.zeros(W->get_size(), m.n_elem);

    return Newmark::initialize();
}

void NonviscousNewmark::assemble_resistance() {
    Newmark::assemble_resistance();

    auto& W = get_domain()->get_factory();

    W->set_sushi(W->get_sushi() + sum(trial_damping, 1) - sum(get_residual(), 1));
}

void NonviscousNewmark::assemble_matrix() {
    Newmark::assemble_matrix();

    auto& W = get_domain()->get_factory();

    const auto damping_diag = C1 * accu(m / (C5 + s));

    for(auto I = 0u; I < W->get_size(); ++I) W->get_stiffness()->unsafe_at(I, I) += damping_diag;
}

int NonviscousNewmark::update_internal(const mat& samurai) {
    trial_damping += C1 * samurai * (m / (C5 + s)).t() - get_residual();

    return SUANPAN_SUCCESS;
}

void NonviscousNewmark::commit_status() {
    current_damping = trial_damping;
    Newmark::commit_status();
}

void NonviscousNewmark::clear_status() {
    trial_damping = current_damping.zeros();
    Newmark::clear_status();
}

void NonviscousNewmark::reset_status() {
    trial_damping = current_damping;
    Newmark::reset_status();
}

void NonviscousNewmark::print() {
    suanpan_info("A NonviscousNewmark solver.\n");
}
