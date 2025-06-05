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
 * @class Newmark
 * @brief A Newmark class defines a solver using Newmark algorithm.
 *
 * `Newmark` algorithm is unconditionally stable if
 * \f{gather}{\alpha\geq\dfrac{1}{4}\left(\dfrac{1}{2}+\beta\right)^2,\qquad\beta\geq\dfrac{1}{2}\f}.
 *
 * There are several choices for solver parameters.
 *
 * Constant Acceleration:
 * \f{gather}{\alpha=\dfrac{1}{4},\qquad\beta=\dfrac{1}{2}\f}.
 *
 * Linear Acceleration:
 * \f{gather}{\alpha=\dfrac{1}{6},\qquad\beta=\dfrac{1}{2}\f}.
 *
 * @author tlc
 * @date 25/08/2017
 * @version 0.1.1
 * @file Newmark.h
 * @addtogroup Integrator
 * @{
 */

#ifndef NEWMARK_H
#define NEWMARK_H

#include "Integrator.h"

class Newmark : public ImplicitIntegrator {
    const double beta;  /**< parameter = .25 */
    const double gamma; /**< parameter = .5 */
protected:
    double C0 = 0., C1 = 0., C2 = 0., C3 = 0., C4 = 0., C5 = 0.; /**< parameters */
public:
    explicit Newmark(unsigned = 0, double = .25, double = .5);

    void assemble_resistance() override;
    void assemble_matrix() override;

    int update_trial_status(bool) override;

    void update_parameter(double) override;

    vec from_incre_velocity(const vec&, const uvec&) override;
    vec from_incre_acceleration(const vec&, const uvec&) override;

    void print() override;
};

#endif

//! @}
