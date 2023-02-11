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
 * @class CustomMises1D
 * @brief A CustomMises1D material class.
 * @author tlc
 * @date 21/01/2023
 * @version 0.1.0
 * @file CustomMises1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CUSTOMMISES1D_H
#define CUSTOMMISES1D_H

#include <Material/Material1D/vonMises/NonlinearMises1D.h>
#include <Toolbox/Expression.h>

class CustomMises1D final : public NonlinearMises1D {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;
    [[nodiscard]] double compute_h(double) const override;
    [[nodiscard]] double compute_dh(double) const override;

    const unsigned k_tag, h_tag;

    ExpressionHolder k_expression, h_expression;

public:
    CustomMises1D(unsigned,   // tag
                  double,     // elastic modulus
                  unsigned,   // isotropic hardening function
                  unsigned,   // kinematic hardening function
                  double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
