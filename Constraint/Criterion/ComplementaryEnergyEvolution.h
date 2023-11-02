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
 * @class ComplementaryEnergyEvolution
 * @brief A ComplementaryEnergyEvolution class.
 *
 * The ComplementaryEnergyEvolution class.
 *
 * @author tlc
 * @date 30/06/2021
 * @version 0.1.0
 * @file ComplementaryEnergyEvolution.h
 * @addtogroup Criterion
 * @{
 */

#ifndef COMPLEMENTARYENERGYEVOLUTION_H
#define COMPLEMENTARYENERGYEVOLUTION_H

#include "EnergyEvolution.h"

class ComplementaryEnergyEvolution final : public EnergyEvolution {
public:
    ComplementaryEnergyEvolution(
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
