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

class BatheTwoStep final : public Integrator {
    enum class FLAG {
        TRAP,
        EULER
    };

    FLAG step_flag = FLAG::TRAP;

    double C0 = 0., C1 = 0., C2 = 0., C3 = 0., C4 = 0., C5 = 0., C6 = 0.;

public:
    using Integrator::Integrator;

    void assemble_resistance() override;
    void assemble_matrix() override;

    int update_trial_status() override;

    void commit_status() override;
    void clear_status() override;

    void update_parameter(double) override;
    void update_compatibility() const override;

    vec from_incre_velocity(const vec&, const uvec&) override;
    vec from_incre_acceleration(const vec&, const uvec&) override;

    void print() override;
};

#endif

//! @}
