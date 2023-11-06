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
 * @class BilinearMises1D
 * @brief A BilinearMises1D material class.
 * @author tlc
 * @date 21/01/2019
 * @version 0.1.0
 * @file BilinearMises1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef BILINEARMISES1D_H
#define BILINEARMISES1D_H

#include <Material/Material1D/vonMises/NonlinearMises1D.h>

struct DataBilinearMises1D {
    const double yield_stress;
    const double isotropic_modulus;
};

class BilinearMises1D final : protected DataBilinearMises1D, public NonlinearMises1D {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;
    [[nodiscard]] double compute_h(double) const override;
    [[nodiscard]] double compute_dh(double) const override;

public:
    BilinearMises1D(
        unsigned,   // tag
        double,     // elastic modulus
        double,     // initial yield stress
        double,     // hardening ratio
        double = 0. // density
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
