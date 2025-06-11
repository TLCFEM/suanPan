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
 * @class BatheExplicit
 * @brief A BatheExplicit class defines a solver using BatheExplicit algorithm.
 *
 * @author tlc
 * @date 03/12/2022
 * @version 0.1.0
 * @file BatheExplicit.h
 * @addtogroup Integrator
 * @{
 */

#ifndef BATHEEXPLICIT_H
#define BATHEEXPLICIT_H

#include "Integrator.h"

class BatheExplicit final : public ExplicitIntegrator {
    enum class FLAG {
        FIRST,
        SECOND
    };

    FLAG step_flag = FLAG::FIRST;

    const double P, Q1, Q2, Q0;
    double DT{0.}, A0{0.}, A1{0.}, A2{0.}, A3{0.}, A4{0.}, A5{0.}, A6{0.}, A7{0.};

protected:
    void update_parameter(double) override;

    int correct_trial_status() override;

public:
    BatheExplicit(unsigned, double);

    [[nodiscard]] bool has_corrector() const override;
    [[nodiscard]] bool time_independent_matrix() const override;

    void update_incre_time(double) override;

    int update_trial_status(bool) override;

    void commit_status() override;
    void clear_status() override;

    void print() override;
};

#endif

//! @}
