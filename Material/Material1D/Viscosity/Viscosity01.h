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
 * @class Viscosity01
 * @brief A 1D Elastic class.
 * @author tlc
 * @date 01/03/2019
 * @version 0.3.0
 * @file Viscosity01.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef VISCOSITY01_H
#define VISCOSITY01_H

#include "NonlinearViscosity.h"

struct DataViscosity01 {
    const double damping;
};

class Viscosity01 final : protected DataViscosity01, public NonlinearViscosity {
    [[nodiscard]] double compute_du(double, double) const override;
    [[nodiscard]] double compute_dv(double, double) const override;
    [[nodiscard]] double compute_damping_coefficient(double, double) const override;

public:
    Viscosity01(unsigned, // tag
                double,   // alpha
                double,   // damp coefficient
                double    // cut-off
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
