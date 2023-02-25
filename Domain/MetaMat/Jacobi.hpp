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
/**
 * @class Jacobi
 * @brief A Jacobi class.
 *
 * @author tlc
 * @date 21/07/2022
 * @version 0.1.0
 * @file Jacobi.hpp
 * @addtogroup Preconditioner
 * @{
 */

#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "Preconditioner.hpp"

template<sp_d data_t> class Jacobi final : public Preconditioner<data_t> {
    const Col<data_t> diag_reciprocal;

public:
    template<typename Container> explicit Jacobi(const Container& in_mat)
        : Preconditioner<data_t>()
        , diag_reciprocal(1. / Col<data_t>(in_mat.diag())) {}

    [[nodiscard]] Col<data_t> apply(const Col<data_t>&) override;
};

template<sp_d data_t> Col<data_t> Jacobi<data_t>::apply(const Col<data_t>& in) {
    Col<data_t> out = in;

    for(auto I = 0llu; I < in.n_elem; I += diag_reciprocal.n_elem) out.subvec(I, arma::size(diag_reciprocal)) %= diag_reciprocal;

    return out;
}

#endif

//! @}
