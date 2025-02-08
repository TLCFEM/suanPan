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
 * @class NLE3D01
 * @brief The NLE3D01 class.
 *
 * @author tlc
 * @date 23/09/2020
 * @version 1.0.0
 * @file NLE3D01.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NLE3D01_H
#define NLE3D01_H

#include "IsotropicNonlinearElastic3D.h"

struct DataNLE3D01 {
    const double bulk; // 9K
    const double ref_strain;
    const double ref_stress;
    const double m;
};

class NLE3D01 final : protected DataNLE3D01, public IsotropicNonlinearElastic3D {
    const double factor_a = (.5 * m + .5) * ref_stress / (1. + m) * pow(ref_strain, -m);
    const double factor_b = (.5 * m - .5) * factor_a;

    vec compute_derivative(double, double) override;

public:
    NLE3D01(
        unsigned, // tag
        double,   // bulk modulus
        double,   // reference strain
        double,   // reference stress
        double,   // m
        double    // density
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
