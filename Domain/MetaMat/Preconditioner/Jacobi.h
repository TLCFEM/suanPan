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
/**
 * @class Jacobi
 * @brief A Jacobi class.
 *
 * @author tlc
 * @date 21/07/2022
 * @version 0.1.0
 * @file Jacobi.h
 * @addtogroup Preconditioner
 * @{
 */

#ifndef JACOBI_H
#define JACOBI_H

#include "Preconditioner.h"

class Jacobi final : public Preconditioner {
    const vec diag_reciprocal;
public:
    template<typename Container> explicit Jacobi(const Container& in_mat)
        : Preconditioner()
        , diag_reciprocal([&] {
            vec t_diag = in_mat.diag();
            return 1. / t_diag.replace(0., t_diag.max());
        }()) {}

    explicit Jacobi(vec&&);

    [[nodiscard]] vec apply(const vec&) const override;

    [[nodiscard]] unique_ptr<Preconditioner> get_copy() const override;
};

#endif

//! @}
