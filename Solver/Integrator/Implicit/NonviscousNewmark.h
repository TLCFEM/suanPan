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
 * @class NonviscousNewmark
 * @brief A NonviscousNewmark class defines a solver using Newmark algorithm.
 *
 * Using exponential function based convolution for global damping model.
 *
 * Reference: 10.1016/j.ymssp.2024.111156
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

class NonviscousNewmark final : public Newmark {
    const cx_vec m, s;

    cx_vec s_para, m_para;

    double accu_para{0.};

    cx_mat current_damping;

protected:
    void update_parameter(double) override;

public:
    NonviscousNewmark(unsigned, double, double, cx_vec&&, cx_vec&&);

    int initialize() override;

    void assemble_resistance() override;
    void assemble_matrix() override;

    void commit_status() override;
    void clear_status() override;

    void print() override;
};

#endif

//! @}
