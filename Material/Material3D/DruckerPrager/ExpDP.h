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
 * @class ExpDP
 * @brief The ExpDP class.
 *
 * Reference:
 *   1. A plastic-damage model for concrete. https://doi.org/10.1016/0020-7683(89)90050-4
 *
 * @author tlc
 * @date 26/10/2018
 * @version 0.1.0
 * @file ExpDP.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef EXPDP_H
#define EXPDP_H

#include "NonlinearDruckerPrager.h"

struct DataExpDP {
    const double cohesion, a, b;
};

class ExpDP final : DataExpDP, public NonlinearDruckerPrager {
    [[nodiscard]] double compute_c(double) const override;
    [[nodiscard]] double compute_dc(double) const override;

public:
    ExpDP(unsigned,   // tag
          double,     // elastic modulus
          double,     // poisson's ratio
          double,     // eta_yield
          double,     // eta_flow
          double,     // xi
          double,     // cohesion
          double,     // a
          double,     // b
          double = 0. // density
        );

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
