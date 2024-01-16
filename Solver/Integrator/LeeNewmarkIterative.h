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
 * @class LeeNewmarkIterative
 * @brief A LeeNewmarkIterative class defines a solver using Newmark algorithm with Lee damping model.
 *
 * Remarks:
 *   1. Only Type 0 is implemented.
 *   2. An iterative algorithm is used. The quadratic convergence is lost.
 *
 * @author tlc
 * @date 15/01/2024
 * @version 0.1.0
 * @file LeeNewmark.h
 * @addtogroup Integrator
 * @{
 */

#ifndef LEENEWMARKITERATIVE_H
#define LEENEWMARKITERATIVE_H

#include "Newmark.h"
#include <Domain/MetaMat/MetaMat.hpp>

class LeeNewmarkIterative final : public Newmark {
    const vec mass_coef, stiffness_coef;

    shared_ptr<MetaMat<double>> current_mass = nullptr;
    shared_ptr<MetaMat<double>> current_stiffness = nullptr;

    void update_damping_force() const;

public:
    LeeNewmarkIterative(unsigned, vec&&, vec&&, double, double);

    [[nodiscard]] int process_constraint() override;
    [[nodiscard]] int process_constraint_resistance() override;

    void assemble_matrix() override;
};

#endif

//! @}
