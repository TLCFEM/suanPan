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
 * @class MultilinearMises1D
 * @brief A MultilinearMises1D material class.
 * @author tlc
 * @date 08/08/2017
 * @version 0.1.0
 * @file MultilinearMises1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef MULTILINEARMISES1D_H
#define MULTILINEARMISES1D_H

#include <Material/Material1D/vonMises/NonlinearMises1D.h>

struct DataMultilinearMises1D {
    const mat backbone;
};

class MultilinearMises1D final : protected DataMultilinearMises1D, public NonlinearMises1D {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;
    [[nodiscard]] double compute_h(double) const override;
    [[nodiscard]] double compute_dh(double) const override;

public:
    MultilinearMises1D(
        unsigned,   // tag
        double,     // elastic modulus
        mat&&,      // backbone
        double = 0. // density
    );

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
