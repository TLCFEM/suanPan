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
 * @class AbsIncreEnergy
 * @brief The AbsIncreEnergy class handles converger test to indicate if the
 * iteration converges.
 *
 * The criterion:
 * \f{gather}{
 * \big|\big|\Delta{}u\cdot{}\Delta{}r\big|\big|_2<E,
 * \f}
 * where \f$E\f$ is the tolerance.
 *
 * @author tlc
 * @date 16/09/2017
 * @version 0.1.0
 * @file AbsIncreEnergy.h
 * @addtogroup Converger
 * @{
 */

#ifndef ABSINCREENERGY_H
#define ABSINCREENERGY_H

#include "Converger.h"

class AbsIncreEnergy final : public Converger {
public:
    explicit AbsIncreEnergy(unsigned = 0, double = 1E-8, unsigned = 7, bool = false);

    unique_ptr<Converger> get_copy() override;

    bool is_converged(unsigned) override;
};

#endif

//! @}
