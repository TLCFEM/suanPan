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
 * @class StrainEnergyEvolution
 * @brief A StrainEnergyEvolution class.
 *
 * The StrainEnergyEvolution class.
 *
 * @author tlc
 * @date 13/09/2020
 * @version 0.1.0
 * @file StrainEnergyEvolution.h
 * @addtogroup Criterion
 * @{
 */

#ifndef STRAINENERGYEVOLUTION_H
#define STRAINENERGYEVOLUTION_H

#include "EnergyEvolution.h"

class StrainEnergyEvolution final : public EnergyEvolution {
public:
    StrainEnergyEvolution(
        unsigned,      // tag
        unsigned,      // step tag
        unsigned,      // incre level
        unsigned,      // final level
        double = 1.,   // centre weight
        unsigned = 2,  // propagation iteration
        unsigned = 10, // reactivation ratio
        double = .5,   // propagation weight
        double = 1E-5  // tolerance
    );

    unique_ptr<Criterion> get_copy() override;
};

#endif

//! @}
