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
 * @class BFGS
 * @brief The (L-)BFGS class defines a solver using BFGS iteration method.
 *
 * The BFGS method is a rank two quasi-Newton method which has a super-linear
 * convergence rate.
 * The (L-)BFGS class supports both conventional BFGS and L-BFGS method
 * which uses limited history information.
 *
 * The Algorithm 7.4 is implemented.
 *
 * References:
 *
 * 1. [Numerical Optimization](https://doi.org/10.1007/978-0-387-40065-5)
 *
 * @author tlc
 * @date 03/05/2025
 * @version 0.4.0
 * @file BFGS.h
 * @addtogroup Solver
 * @{
 */

#ifndef BFGS_H
#define BFGS_H

#include <Solver/Solver.h>
#include <deque>

class BFGS final : public Solver {
    std::deque<vec> hist_s, hist_y;
    std::deque<double> hist_rho;
    std::vector<double> hist_alpha;

    const unsigned max_hist;

public:
    explicit BFGS(unsigned = 0, unsigned = 20);

    int analyze() override;

    void print() override;
};

#endif

//! @}
