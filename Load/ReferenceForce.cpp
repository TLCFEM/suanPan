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

#include "ReferenceForce.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

ReferenceForce::ReferenceForce(const unsigned T, const unsigned S, const double L, uvec&& N, const unsigned D)
    : Load(T, S, 0, std::forward<uvec>(N), uvec{D}, L) {}

int ReferenceForce::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();

    reference_load.zeros(W->get_size());

    for(const auto I : get_nodal_active_dof(D)) reference_load(I) = pattern;

    return SUANPAN_SUCCESS;
}
