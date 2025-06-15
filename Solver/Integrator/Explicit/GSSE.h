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
 * @class GSSE
 * @brief General explicit integration.
 *
 * References:
 *     1. [10.1002/nme.6574](https://doi.org/10.1002/nme.6574)
 *
 * @author tlc
 * @date 15/06/2025
 * @version 0.1.0
 * @file GSSE.h
 * @addtogroup Integrator
 * @{
 */

#ifndef GSSE_H
#define GSSE_H

#include "../Integrator.h"

class GSSE : public ExplicitIntegrator {
    double DT{0.};

    const double C;
    const double UA0, VA0;
    const double UB, VB, AB;

protected:
    void update_parameter(double) override;

    [[nodiscard]] int process_load_impl(bool) override;
    [[nodiscard]] int process_constraint_impl(bool) override;

    [[nodiscard]] bool has_corrector() const override;

    int correct_trial_status() override;

public:
    GSSE(unsigned, double, double);
    GSSE(unsigned, double);

    int update_trial_status(bool) override;

    vec from_incre_acceleration(const vec&, const uvec&) override;
    vec from_total_acceleration(const vec&, const uvec&) override;

    void print() override;
};

class ICL final : public GSSE {
public:
    ICL(unsigned, double);
};

#endif

//! @}
