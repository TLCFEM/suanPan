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
 * @class Newton
 * @brief A Newton class defines a solver using Newton--Raphson iteration.
 * @author tlc
 * @date 27/08/2017
 * @version 0.1.2
 * @file Newton.h
 * @addtogroup Solver
 * @{
 */

#ifndef NEWTON_H
#define NEWTON_H

#include <Solver/Solver.h>

class Newton final : public Solver {
    const bool initial_stiffness;

public:
    explicit Newton(unsigned = 0, bool = false);

    int analyze() override;

    void print() override;
};

#endif

//! @}
