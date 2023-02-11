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
 * @class CustomDegradation
 * @brief The CustomDegradation class.
 *
 * @author tlc
 * @date 24/01/2023
 * @version 0.1.0
 * @file CustomDegradation.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CUSTOMDEGRADATION_H
#define CUSTOMDEGRADATION_H

#include "Degradation.h"
#include <Toolbox/Expression.h>

class CustomDegradation final : public Degradation {
    const unsigned expression_tag;

    ExpressionHolder expression;

    [[nodiscard]] podarray<double> compute_degradation(double) const override;

public:
    CustomDegradation(unsigned, // unique tag
                      unsigned, // material tag
                      unsigned  // expression tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
