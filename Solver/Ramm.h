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
 * @class Ramm
 * @brief A Ramm class defines a solver using Ramm's version of arc--length method.
 * @author tlc
 * @date 27/07/2017
 * @version 0.1.0
 * @file Ramm.h
 * @addtogroup Solver
 * @{
 */

#ifndef RAMM_H
#define RAMM_H

#include <Solver/Solver.h>

class Ramm final : public Solver {
    double arc_length = 1.;

public:
    using Solver::Solver;

    int analyze() override;

    void set_step_size(double) override;

    void print() override;
};

#endif

//! @}
