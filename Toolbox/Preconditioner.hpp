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

#ifndef PRECONDITIONER_HPP
#define PRECONDITIONER_HPP

#include <suanPan.h>

template<typename T> concept IsPreconditioner = requires(T t, const vec& x) { t.apply(x); };

template<sp_d data_t> class DiagonalPreconditioner {
    vec diag_reciprocal;
public:
    explicit DiagonalPreconditioner(const Mat<data_t>& in_mat)
        : diag_reciprocal(1. / in_mat.diag()) {}

    [[nodiscard]] Col<data_t> apply(const Col<data_t>& in) const { return diag_reciprocal % in; }
};

#endif
