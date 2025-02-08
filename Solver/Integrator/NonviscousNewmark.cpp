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

#include "NonviscousNewmark.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

NonviscousNewmark::NonviscousNewmark(const unsigned T, const double A, const double B, cx_vec&& M, cx_vec&& S)
    : Newmark(T, A, B)
    , m(std::move(M))
    , s(std::move(S)) {}

int NonviscousNewmark::initialize() {
    auto& W = get_domain()->get_factory();

    current_damping.zeros(W->get_size(), m.n_elem);

    return Newmark::initialize();
}

void NonviscousNewmark::assemble_resistance() {
    Newmark::assemble_resistance();

    auto& W = get_domain()->get_factory();

    const vec trial_damping_force = real(current_damping * s_para + accu_para * (W->get_current_velocity() + W->get_trial_velocity()));

    W->update_trial_damping_force_by(trial_damping_force);

    W->update_sushi_by(trial_damping_force);
}

void NonviscousNewmark::assemble_matrix() {
    Newmark::assemble_matrix();

    auto& W = get_domain()->get_factory();
    const auto& t_stiffness = W->get_stiffness();

    const auto damping_diag = C1 * accu_para;

    if(const auto t_scheme = W->get_storage_scheme(); StorageScheme::SPARSE == t_scheme || StorageScheme::SPARSESYMM == t_scheme) for(auto I = 0u; I < W->get_size(); ++I) t_stiffness->unsafe_at(I, I) += damping_diag;
    else suanpan::for_each(W->get_size(), [&](const unsigned I) { t_stiffness->unsafe_at(I, I) += damping_diag; });
}

void NonviscousNewmark::update_parameter(const double DT) {
    Newmark::update_parameter(DT);

    const cx_vec t_para = 2. / DT + s;
    s_para = (t_para - 2. * s) / t_para;
    m_para = m / t_para;
    accu_para = accu(m_para).real();
}

void NonviscousNewmark::commit_status() {
    auto& W = get_domain()->get_factory();

    current_damping *= diagmat(s_para);
    current_damping += (W->get_current_velocity() + W->get_trial_velocity()) * m_para.t();

    Newmark::commit_status();
}

void NonviscousNewmark::clear_status() {
    current_damping.zeros();
    Newmark::clear_status();
}

void NonviscousNewmark::print() {
    suanpan_info("A NonviscousNewmark solver.\n");
}
