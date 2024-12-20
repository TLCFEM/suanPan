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

#include "GroupNodalForce.h"
#include <Domain/DomainBase.h>

GroupNodalForce::GroupNodalForce(const unsigned T, const unsigned S, const double L, uvec&& N, const unsigned D, const unsigned AT)
    : GroupLoad(std::move(N))
    , NodalForce(T, S, L, uvec{}, uvec{D}, AT) {}

GroupNodalForce::GroupNodalForce(const unsigned T, const unsigned S, const double L, uvec&& N, uvec&& D, const unsigned AT)
    : GroupLoad(std::move(N))
    , NodalForce(T, S, L, uvec{}, std::move(D), AT) {}

int GroupNodalForce::initialize(const shared_ptr<DomainBase>& D) {
    node_encoding = update_object_tag(D);

    return NodalForce::initialize(D);
}
