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
 * @class NonviscousNewmark
 * @brief A NonviscousNewmark class defines a solver using Newmark algorithm.
 *
 * Using exponential function based convolution for global damping model.
 *
 * @author tlc
 * @date 18/03/2023
 * @version 0.1.0
 * @file NonviscousNewmark.h
 * @addtogroup Integrator
 * @{
 */

#ifndef NONVISCOUSNEWMARK_H
#define NONVISCOUSNEWMARK_H

#include "Newmark.h"

class NonviscousNewmark : public Newmark {
    const vec m, s;

    mat trial_damping, current_damping;

    mat get_residual() const;

public:
    explicit NonviscousNewmark(unsigned, double, double, vec&&, vec&&);

    int initialize() override;

    void assemble_resistance() override;
    void assemble_matrix() override;

    int update_internal(const mat&) override;

    void commit_status() override;
    void clear_status() override;
    void reset_status() override;

    void print() override;
};

#endif

//! @}
