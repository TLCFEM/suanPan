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
 * @class ExpJ2
 * @brief The ExpJ2 class.
 *
 * The ExpJ2 class defines a 3D material model using the J2 yield criterion.
 * There is no kinematic hardening defined so parameter \f$beta\f$ only controls
 * the amount of isotropic hardening, which is defined using the following expression
 * \f{align}{
 * \sigma=\sigma_y\left(\left(1+a\right)\exp\left(-b\varepsilon_p\right)-a\exp\left(-2b\varepsilon_p\right)\right)
 * \f}
 * where \f$a\f$ and \f$b\f$ are two parameters.
 *
 * If \f$a>1\f$, hardening are defined, the ratio \f$\sigma_m/\sigma_y=\dfrac{(1+a)^2}{4a}\f$.
 *
 * Reference:
 * 1. A plastic-damage model for concrete. https://doi.org/10.1016/0020-7683(89)90050-4
 *
 * @author tlc
 * @date 26/10/2018
 * @version 0.1.0
 * @file ExpJ2.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef EXPJ2_H
#define EXPJ2_H

#include "NonlinearJ2.h"

struct DataExpJ2 {
    const double yield_stress;
    const double a, b;
};

class ExpJ2 final : DataExpJ2, public NonlinearJ2 {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;
    [[nodiscard]] double compute_h(double) const override;
    [[nodiscard]] double compute_dh(double) const override;

public:
    ExpJ2(unsigned,   // tag
          double,     // elastic modulus
          double,     // poisson's ratio
          double,     // yield stress
          double,     // a
          double,     // b
          double = 0. // density
    );

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
