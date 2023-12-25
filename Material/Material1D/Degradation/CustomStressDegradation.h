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
 * @class CustomStressDegradation
 * @brief The CustomStressDegradation class.
 *
 * @author tlc
 * @date 02/09/2023
 * @version 0.1.0
 * @file CustomStressDegradation.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CUSTOMSTRESSDEGRADATION_H
#define CUSTOMSTRESSDEGRADATION_H

#include "Degradation.h"
#include <Toolbox/Expression.h>
#include <Toolbox/ResourceHolder.h>

class CustomStressDegradation final : public StressDegradation {
    const unsigned positive_expression_tag, negative_expression_tag;

    ResourceHolder<Expression> positive_expression, negative_expression;

    [[nodiscard]] vec compute_positive_degradation(double) const override;
    [[nodiscard]] vec compute_negative_degradation(double) const override;

public:
    CustomStressDegradation(
        unsigned, // unique tag
        unsigned, // material tag
        unsigned, // expression tag
        unsigned  // expression tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
