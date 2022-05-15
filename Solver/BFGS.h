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
 * @class BFGS
 * @brief The (L-)BFGS class defines a solver using BFGS iteration method.
 *
 * The BFGS method is a rank two quasi--Newton method which has a super linear
 * convergence rate. The (L-)BFGS class supports both conventional BFGS and
 * L-BFGS method which uses limited history information.
 *
 * \f{gather}{
 * K_{n+1}^{-1}=\left(I-\dfrac{\Delta{}UR^T}{R^T\Delta{}U}\right)K_n^{-1}\left(I-\dfrac{R\Delta{}U^T}{R^T\Delta{}U}\right)+\dfrac{\Delta{}U\Delta{}U^T}{R^T\Delta{}U}.
 * \f}
 *
 * The \f$I\f$ is identity matrix. The \f$\Delta{}U\f$ is current displacement
 * increment. The \f$R\f$ is current residual. For brevity, in both terms, the
 * subscript \f$n\f$ representing current step is dropped.
 *
 * @author tlc
 * @date 06/07/2018
 * @version 0.3.0
 * @file BFGS.h
 * @addtogroup Solver
 * @{
 */

#ifndef BFGS_H
#define BFGS_H

#include <Solver/Solver.h>
#include <deque>

using std::deque;
using std::vector;

class BFGS final : public Solver {
    deque<vec> hist_ninja, hist_residual;
    deque<double> hist_factor;
    vector<double> alpha;

    const unsigned max_hist;

public:
    explicit BFGS(unsigned = 0, unsigned = 30);

    int analyze() override;

    void print() override;
};

#endif

//! @}
