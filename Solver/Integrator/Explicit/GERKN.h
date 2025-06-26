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
 * @class GERKN
 * @brief Generalized Explicit Runge-Kutta-Nystrom time integration.
 *
 * References:
 *     1. [10.1002/nme.7658](https://doi.org/10.1002/nme.7658)
 *
 * @author tlc
 * @date 08/06/2025
 * @version 0.1.0
 * @file GERKN.h
 * @addtogroup Integrator
 * @{
 */

#ifndef GERKN_H
#define GERKN_H

#include "../Integrator.h"

class GERKN : public ExplicitIntegrator {
    enum class FLAG {
        FIRST,
        SECOND
    };

    FLAG step_flag = FLAG::FIRST;

    double DT{0.};

protected:
    void update_parameter(double) override;

    [[nodiscard]] int process_load_impl(bool) override;
    [[nodiscard]] int process_constraint_impl(bool) override;

    [[nodiscard]] bool has_corrector() const override;

    int correct_trial_status() override;

    double C1{0.}, C2{0.};

    double UA10{0.}, UA20{0.}, UA21{0.};
    double UB0{0.}, UB1{0.}, UB2{0.};

    double VA10{0.}, VA20{0.}, VA21{0.};
    double VB0{0.}, VB1{0.}, VB2{0.};

    double AB0{0.}, AB1{0.}, AB2{0.};

public:
    using ExplicitIntegrator::ExplicitIntegrator;

    void update_incre_time(double) override;

    int update_trial_status(bool) override;

    void commit_status() override;
    void clear_status() override;

    vec from_incre_acceleration(const vec&, const uvec&) override;
    vec from_total_acceleration(const vec&, const uvec&) override;
};

class WAT2 final : public GERKN {
public:
    WAT2(unsigned, double);

    void print() override;
};

#endif

//! @}
