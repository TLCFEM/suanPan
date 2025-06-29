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
 * @class CustomHoffman
 * @brief The CustomHoffman class.
 *
 * @author tlc
 * @date 16/01/2023
 * @version 0.2.0
 * @file CustomHoffman.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef CUSTOMHOFFMAN_H
#define CUSTOMHOFFMAN_H

#include "NonlinearOrthotropic.h"

#include <Toolbox/Expression.h>
#include <Toolbox/ResourceHolder.h>

class CustomHoffman final : public NonlinearOrthotropic {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;

    const unsigned k_tag;

    ResourceHolder<Expression> k_expression;

public:
    CustomHoffman(
        unsigned,   // tag
        vec&&,      // elastic modulus
        vec&&,      // poissons ratio
        vec&&,      // sigma
        unsigned,   // hardening function tag
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
