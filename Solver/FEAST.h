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
 * @class FEAST
 * @brief A FEAST class defines a solver using FEAST iteration.
 * @author tlc
 * @date 12/04/2021
 * @version 0.1.0
 * @file FEAST.h
 * @addtogroup Solver
 * @{
 */

#ifndef FEAST_H
#define FEAST_H

#include <Solver/Solver.h>

template<sp_d T> class Factory;
using LongFactory = Factory<double>;

class FEAST final : public Solver {
    static char UPLO;

    const bool quadratic = false;

    const unsigned eigen_num;
    const double centre, radius;

    [[nodiscard]] int linear_solve(const shared_ptr<LongFactory>&) const;
    [[nodiscard]] int quadratic_solve(const shared_ptr<LongFactory>&) const;

public:
    FEAST(unsigned, unsigned, double, double, bool);

    int initialize() override;

    int analyze() override;

    void print() override;
};

#endif

//! @}
