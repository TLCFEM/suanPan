/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
 * @class ExpGurson1D
 * @brief The ExpGurson1D class.
 *
 * @author tlc
 * @date 06/01/2020
 * @version 0.1.0
 * @file ExpGurson1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef EXPGURSON1D_H
#define EXPGURSON1D_H

#include "NonlinearGurson1D.h"

struct DataExpGurson1D {
    const double yield_stress, n;
};

class ExpGurson1D final : DataExpGurson1D, public NonlinearGurson1D {
    static const unsigned max_iteration;

    const double para_c = 3. * elastic_modulus / (2. + 2. * poissons_ratio) / yield_stress;

    [[nodiscard]] vec compute_hardening(double) const override;

public:
    ExpGurson1D(unsigned,   // tag
                double,     // elastic modulus
                double,     // poisson's ratio
                double,     // yield stress
                double,     // n
                double,     // q1
                double,     // q2
                double,     // fn
                double,     // sn
                double,     // en
                double = 0. // density
    );

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
