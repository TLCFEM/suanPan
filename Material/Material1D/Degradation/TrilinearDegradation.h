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
 * @class TrilinearDegradation
 * @brief The TrilinearDegradation class.
 *
 * @author tlc
 * @date 11/05/2019
 * @version 0.1.0
 * @file TrilinearDegradation.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef TRILINEARDEGRADATION_H
#define TRILINEARDEGRADATION_H

#include "Degradation.h"

struct DataTrilinearDegradation {
    const double s_strain;
    const double e_strain;
    const double e_damage;
    const double slope = (1. - e_damage) / (s_strain - e_strain);
};

class TrilinearDegradation final : DataTrilinearDegradation, public Degradation {
    [[nodiscard]] podarray<double> compute_degradation(double) const override;

public:
    TrilinearDegradation(unsigned, // unique tag
                         unsigned, // material tag
                         double,   // start strain
                         double,   // end strain
                         double    // end level
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
