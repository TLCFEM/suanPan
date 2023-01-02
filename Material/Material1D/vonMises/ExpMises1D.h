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
 * @class ExpMises1D
 * @brief A ExpMises1D material class.
 * @author tlc
 * @date 21/01/2019
 * @version 0.1.0
 * @file ExpMises1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef EXPMISES1D_H
#define EXPMISES1D_H

#include <Material/Material1D/vonMises/NonlinearMises1D.h>

struct DataExpMises1D {
    const double yield_stress;
    const double a, b, c;
};

class ExpMises1D final : DataExpMises1D, public NonlinearMises1D {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;
    [[nodiscard]] double compute_h(double) const override;
    [[nodiscard]] double compute_dh(double) const override;

public:
    explicit ExpMises1D(unsigned,   // tag
                        double,     // elastic modulus
                        double,     // initial yield stress
                        double,     // a
                        double,     // b
                        double,     // c
                        double = 0. // density
        );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
