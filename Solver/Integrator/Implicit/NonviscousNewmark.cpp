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

#include "NonviscousNewmark.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

vec NonviscousNewmark::target_field() const {
    auto& W = get_domain()->get_factory();

    return W->get_current_velocity() + W->get_trial_velocity();
}

void NonviscousNewmark::assemble_effective_matrix() {
    auto& W = get_domain()->get_factory();
    const auto& t_stiffness = W->get_stiffness();

    if(W->is_nlgeom()) t_stiffness += W->get_geometry();

    t_stiffness += C0 * W->get_mass();

    t_stiffness += W->is_nonviscous() ? C1 * (W->get_damping() + W->get_nonviscous()) : C1 * W->get_damping();

    const auto damping_diag = C1 * accu_para;

    if(const auto t_scheme = W->get_storage_scheme(); StorageScheme::SPARSE == t_scheme || StorageScheme::SPARSESYMM == t_scheme)
        for(auto I = 0u; I < W->get_size(); ++I) t_stiffness->unsafe_at(I, I) += damping_diag;
    else suanpan::for_each(W->get_size(), [&](const unsigned I) { t_stiffness->unsafe_at(I, I) += damping_diag; });
}

void NonviscousNewmark::print() {
    suanpan_info("A NonviscousNewmark solver using convolution based nonviscous damping. doi:10.1016/j.ymssp.2024.111156\n");
}
