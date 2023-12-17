/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class CustomJ2
 * @brief The CustomJ2 class.
 *
 * @author tlc
 * @date 16/01/2023
 * @version 0.1.0
 * @file CustomJ2.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef CUSTOMJ2_H
#define CUSTOMJ2_H

#include "NonlinearJ2.h"
#include <Toolbox/Expression.h>
#include <Toolbox/ResourceHolder.h>

class CustomJ2 final : public NonlinearJ2 {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;
    [[nodiscard]] double compute_h(double) const override;
    [[nodiscard]] double compute_dh(double) const override;

    const unsigned k_tag, h_tag;

    ResourceHolder<Expression> k_expression, h_expression;

public:
    CustomJ2(
        unsigned,   // tag
        double,     // elastic modulus
        double,     // poisson's ratio
        unsigned,   // isotropic hardening function
        unsigned,   // kinematic hardening function
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
