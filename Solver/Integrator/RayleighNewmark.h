/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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
 * @class RayleighNewmark
 * @brief A RayleighNewmark class defines a solver using Newmark algorithm with Rayleigh damping model.
 *
 * `RayleighNewmark` algorithm is unconditionally stable if
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
 * @file RayleighNewmark.h
 * @addtogroup Integrator
 * @{
 */

#ifndef RAYLEIGHNEWMARK_H
#define RAYLEIGHNEWMARK_H

#include <Solver/Integrator/Newmark.h>

class RayleighNewmark final : public Newmark {
	const double damping_alpha;
	const double damping_beta;
	const double damping_zeta;
public:
	explicit RayleighNewmark(unsigned, double, double, double, double = .25, double = .5);

	void assemble_resistance() override;
};

#endif

//! @}
