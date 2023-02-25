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
 * @class CustomViscosity
 * @brief A 1D Viscosity class.
 * @author tlc
 * @date 23/01/2023
 * @version 0.1.0
 * @file CustomViscosity.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CUSTOMVISCOSITY_H
#define CUSTOMVISCOSITY_H

#include "NonlinearViscosity.h"
#include <Toolbox/Expression.h>
#include <Toolbox/ResourceHolder.h>

class CustomViscosity final : public NonlinearViscosity {
    [[nodiscard]] double compute_du(double, double) const override;
    [[nodiscard]] double compute_dv(double, double) const override;
    [[nodiscard]] double compute_damping_coefficient(double, double) const override;

    const unsigned expression_tag;

    ResourceHolder<Expression> expression;

public:
    CustomViscosity(unsigned, // tag
                    unsigned  // expression tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
