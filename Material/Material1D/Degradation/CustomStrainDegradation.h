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
 * @class CustomStrainDegradation
 * @brief The CustomStrainDegradation class.
 *
 * @author tlc
 * @date 24/01/2023
 * @version 0.1.0
 * @file CustomStrainDegradation.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CUSTOMSTRAINDEGRADATION_H
#define CUSTOMSTRAINDEGRADATION_H

#include "Degradation.h"
#include <Toolbox/Expression.h>
#include <Toolbox/ResourceHolder.h>

class CustomStrainDegradation final : public StrainDegradation {
    const unsigned positive_expression_tag, negative_expression_tag;

    ResourceHolder<Expression> positive_expression, negative_expression;

    [[nodiscard]] vec compute_positive_degradation(double) const override;
    [[nodiscard]] vec compute_negative_degradation(double) const override;

public:
    CustomStrainDegradation(
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
