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
 * @class BilinearViscosity
 * @brief A 1D Viscosity class.
 * @author tlc
 * @date 01/09/2020
 * @version 0.1.0
 * @file BilinearViscosity.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef BILINEARVISCOSITY_H
#define BILINEARVISCOSITY_H

#include "NonlinearViscosity.h"

struct DataBilinearViscosity {
    const double damping;
    const double yield_stress;
    const double hardening;
};

class BilinearViscosity final : protected DataBilinearViscosity, public NonlinearViscosity {
    const double yield_strain = yield_stress / damping;

    [[nodiscard]] double compute_du(double, double) const override;
    [[nodiscard]] double compute_dv(double, double) const override;
    [[nodiscard]] double compute_damping_coefficient(double, double) const override;

public:
    BilinearViscosity(
        unsigned, // tag
        double,   // damping
        double,   // yield stress
        double    // hardening ratio
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
