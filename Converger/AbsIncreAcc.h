/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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
 * @class AbsIncreAcc
 * @brief The AbsIncreAcc class that handles converger test to indicate if the
 * iteration converges.
 * @author tlc
 * @date 03/12/2022
 * @version 0.1.0
 * @file AbsIncreAcc.h
 * @addtogroup Converger
 * @{
 */

#ifndef ABSINCREACC_H
#define ABSINCREACC_H

#include "Converger.h"

class AbsIncreAcc final : public Converger {
public:
    explicit AbsIncreAcc(unsigned = 0, double = 1E-8, unsigned = 7, bool = false);

    unique_ptr<Converger> get_copy() override;

    bool is_converged(unsigned) override;
};

#endif

//! @}
