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
 * @class OALTS
 * @brief A OALTS class defines a solver using OALTS algorithm.
 *
 * doi:10.1002/nme.6188
 *
 * @author tlc
 * @date 05/12/2022
 * @version 0.1.0
 * @file OALTS.h
 * @addtogroup Integrator
 * @{
 */

#ifndef OALTS_H
#define OALTS_H

#include "../Integrator.h"

class OALTS final : public ImplicitIntegrator {
    const double A1, A2, B0, B1, B2, B10, B20;

    double DT{0.}, P1{0.}, P2{0.}, P3{0.};

    bool if_starting = true;

protected:
    void update_parameter(double) override;

public:
    OALTS(unsigned, double);

    [[nodiscard]] bool time_independent_matrix() const override;

    void assemble_resistance() override;
    void assemble_matrix() override;

    int update_trial_status(bool) override;

    void commit_status() override;
    void clear_status() override;

    vec from_incre_velocity(const vec&, const uvec&) override;
    vec from_incre_acceleration(const vec&, const uvec&) override;
    vec from_total_velocity(const vec&, const uvec&) override;
    vec from_total_acceleration(const vec&, const uvec&) override;

    void print() override;
};

#endif

//! @}
