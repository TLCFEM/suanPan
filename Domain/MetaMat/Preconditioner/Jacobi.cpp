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

#include "Jacobi.h"

Jacobi::Jacobi(vec&& in_diag)
    : Preconditioner()
    , diag_reciprocal(1. / in_diag.replace(0., in_diag.max())) {}

vec Jacobi::apply(const vec& in) const {
    vec out = in;

    for(auto I = 0llu; I < in.n_elem; I += diag_reciprocal.n_elem) out.subvec(I, size(diag_reciprocal)) %= diag_reciprocal;

    return out;
}

unique_ptr<Preconditioner> Jacobi::get_copy() const { return make_unique<Jacobi>(*this); }
