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
 * @class Dhakal
 * @brief The Dhakal class.
 *
 * @author tlc
 * @date 14/07/2019
 * @version 0.1.0
 * @file Dhakal.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef DHAKAL_H
#define DHAKAL_H

#include "Degradation.h"

struct DataDhakal {
    const double yield_strain, inter_strain, inter_factor;
    const double slope = (inter_factor - 1.) / (inter_strain - yield_strain);
    const double final_strain = (inter_factor - .2) * 50. * yield_strain + inter_strain;
};

class Dhakal final : DataDhakal, public Degradation {
    [[nodiscard]] podarray<double> compute_degradation(double) const override;

public:
    Dhakal(unsigned, // unique tag
           unsigned, // material tag
           double,   // yield strain
           double    // parameter
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
