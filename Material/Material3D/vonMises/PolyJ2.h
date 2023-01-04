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
 * @class PolyJ2
 * @brief The PolyJ2 class.
 *
 * @author tlc
 * @date 25/07/2018
 * @version 0.1.0
 * @file PolyJ2.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef POLYJ2_H
#define POLYJ2_H

#include "NonlinearJ2.h"

struct DataPolyJ2 {
    const vec poly_para;
};

class PolyJ2 final : DataPolyJ2, public NonlinearJ2 {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;
    [[nodiscard]] double compute_h(double) const override;
    [[nodiscard]] double compute_dh(double) const override;

public:
    PolyJ2(unsigned,   // tag
           double,     // elastic modulus
           double,     // poisson's ratio
           vec&&,      // poly parameter
           double = 0. // density
    );

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
