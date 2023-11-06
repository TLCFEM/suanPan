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
 * @class Arnoldi
 * @brief A Arnoldi class defines a solver using Arnoldi iteration.
 * @author tlc
 * @date 19/10/2017
 * @version 0.1.0
 * @file Arnoldi.h
 * @addtogroup Solver
 * @{
 */

#ifndef ARNOLDI_H
#define ARNOLDI_H

#include <Solver/Solver.h>

class Arnoldi final : public Solver {
    const unsigned eigen_num;
    const char eigen_type;

public:
    explicit Arnoldi(
        unsigned = 0, // unique solver tag
        unsigned = 1, // number of eigenvalues
        char = 'S'    // type
    );

    int initialize() override;

    int analyze() override;

    void print() override;
};

#endif

//! @}
