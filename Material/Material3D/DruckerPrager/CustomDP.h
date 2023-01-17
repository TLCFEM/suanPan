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
 * @class CustomDP
 * @brief The CustomDP class.
 *
 * @author tlc
 * @date 21/01/2019
 * @version 0.1.0
 * @file CustomDP.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef CUSTOMDP_H
#define CUSTOMDP_H

#include "NonlinearDruckerPrager.h"
#include <Toolbox/Expression.h>

class CustomDP final : public NonlinearDruckerPrager {
    [[nodiscard]] double compute_c(double) const override;
    [[nodiscard]] double compute_dc(double) const override;

    const unsigned c_tag;

    shared_ptr<Expression> c_expression;

public:
    CustomDP(unsigned,   // tag
             double,     // elastic modulus
             double,     // poisson's ratio
             double,     // eta_yield (hydrostatic stress related)
             double,     // eta_flow (dilatancy angle related)
             double,     // xi (cohesion related)
             unsigned,   // expression tag
             double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
