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
 * @class TrilinearStrainDegradation
 * @brief The TrilinearStrainDegradation class.
 *
 * @author tlc
 * @date 11/05/2019
 * @version 0.1.0
 * @file TrilinearStrainDegradation.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef TRILINEARSTRAINDEGRADATION_H
#define TRILINEARSTRAINDEGRADATION_H

#include "Degradation.h"

struct DataTrilinearStrainDegradation {
    const double s_strain;
    const double e_strain;
    const double e_damage;
    const double slope = (1. - e_damage) / (s_strain - e_strain);
};

class TrilinearStrainDegradation final : protected DataTrilinearStrainDegradation, public StrainDegradation {
    [[nodiscard]] vec compute_positive_degradation(double) const override;
    [[nodiscard]] vec compute_negative_degradation(double) const override;

public:
    TrilinearStrainDegradation(unsigned, // unique tag
                               unsigned, // material tag
                               double,   // start strain
                               double,   // end strain
                               double    // end level
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
