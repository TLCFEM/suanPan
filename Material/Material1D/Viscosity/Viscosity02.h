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
 * @class Viscosity02
 * @brief A 1D Viscosity class.
 * @author tlc
 * @date 01/03/2019
 * @version 0.2.0
 * @file Viscosity02.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef VISCOSITY02_H
#define VISCOSITY02_H

#include "NonlinearViscosity.h"

struct DataViscosity02 {
    const double damping_a, damping_b, damping_c, damping_d;
    const double gap_a, gap_b;
};

class Viscosity02 final : DataViscosity02, public NonlinearViscosity {
    [[nodiscard]] double compute_du(double, double) const override;
    [[nodiscard]] double compute_dv(double, double) const override;
    [[nodiscard]] double compute_damping_coefficient(double, double) const override;

public:
    Viscosity02(unsigned, // tag
                double,   // alpha
                double,   // damp_a
                double,   // damp_b
                double,   // damp_c
                double,   // damp_d
                double,   // gap_a
                double,   // gap_b
                double    // cut-off
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
