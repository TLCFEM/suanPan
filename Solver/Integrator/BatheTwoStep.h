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
 * @class BatheTwoStep
 * @brief A BatheTwoStep class defines a solver using BatheTwoStep algorithm.
 *
 * @author tlc
 * @date 13/05/2020
 * @version 0.1.0
 * @file BatheTwoStep.h
 * @addtogroup Integrator
 * @{
 */

#ifndef BATHETWOSTEP_H
#define BATHETWOSTEP_H

#include "Integrator.h"

class BatheTwoStep final : public ImplicitIntegrator {
    enum class FLAG {
        TRAP,
        EULER
    };

    FLAG step_flag = FLAG::TRAP;

    const double GM;

    const double Q1, Q2, Q0, Q02 = Q0 / Q2, Q12 = Q1 / Q2;

    double P0{0.}, P1{0.}, P2{0.}, P3{0.}, P4{0.}, P5{0.}, P6{0.}, P7{0.}, P8{0.}, P9{0.};

public:
    BatheTwoStep(unsigned, double, double);

    void assemble_resistance() override;
    void assemble_matrix() override;

    void update_incre_time(double) override;

    int update_trial_status() override;

    void commit_status() override;
    void clear_status() override;

    void update_parameter(double) override;

    vec from_incre_velocity(const vec&, const uvec&) override;
    vec from_incre_acceleration(const vec&, const uvec&) override;
    vec from_total_velocity(const vec&, const uvec&) override;
    vec from_total_acceleration(const vec&, const uvec&) override;

    void print() override;
};

#endif

//! @}
