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
 * @class ExpHoffman
 * @brief The ExpHoffman class.
 *
 * @author tlc
 * @date 20/02/2019
 * @version 0.1.0
 * @file ExpHoffman.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef EXPHOFFMAN_H
#define EXPHOFFMAN_H

#include "NonlinearHoffman.h"

struct DataExpHoffman {
    const double a, b;
};

class ExpHoffman final : DataExpHoffman, public NonlinearHoffman {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;

public:
    ExpHoffman(unsigned,   // tag
               vec&&,      // elastic modulus
               vec&&,      // poissons ratio
               vec&&,      // sigma
               double,     // a
               double,     // b
               double = 0. // density
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
