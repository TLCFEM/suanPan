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
/**
 * @class BilinearHoffman
 * @brief The BilinearHoffman class.
 *
 * @author tlc
 * @date 20/01/2019
 * @version 0.2.0
 * @file BilinearHoffman.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef BILINEARHOFFMAN_H
#define BILINEARHOFFMAN_H

#include "NonlinearHoffman.h"

struct DataBilinearHoffman {
    const double hardening_modulus;
};

class BilinearHoffman : protected DataBilinearHoffman, public NonlinearHoffman {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;

public:
    BilinearHoffman(
        unsigned,   // tag
        vec&&,      // elastic modulus
        vec&&,      // poissons ratio
        vec&&,      // sigma
        double,     // hardening ratio
        double = 0. // density
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
