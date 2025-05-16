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
 * @class WilsonPenzienNewmark
 * @brief A WilsonPenzienNewmark class defines a solver using Newmark algorithm with Wilson-Penzien damping model.
 * @author tlc
 * @date 05/06/2020
 * @version 0.1.1
 * @file WilsonPenzienNewmark.h
 * @addtogroup Integrator
 * @{
 */

#ifndef WILSONPENZIENNEWMARK_H
#define WILSONPENZIENNEWMARK_H

#include <Solver/Integrator/Newmark.h>

class WilsonPenzienNewmark final : public Newmark {
    bool first_iteration = true;

    const vec damping_ratio;

    mat theta;
    vec beta;

public:
    WilsonPenzienNewmark(unsigned, vec&&, double = .25, double = .5);

    int initialize() override;

    [[nodiscard]] int process_constraint() override;

    int solve(mat&, const mat&) override;
    int solve(mat&, const sp_mat&) override;
    int solve(mat&, mat&&) override;
    int solve(mat&, sp_mat&&) override;

    void commit_status() override;
    void clear_status() override;
    void reset_status() override;

    void assemble_resistance() override;
    void assemble_matrix() override;

    void print() override;
};

#endif

//! @}
