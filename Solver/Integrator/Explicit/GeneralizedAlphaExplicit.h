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
 * @class GeneralizedAlphaExplicit
 * @brief A GeneralizedAlphaExplicit class defines a solver using GeneralizedAlphaExplicit algorithm.
 *
 * References:
 *     1. [10.1002/nme.6574](https://doi.org/10.1002/nme.6574)
 *
 * @author tlc
 * @date 03/12/2022
 * @version 0.1.0
 * @file GeneralizedAlphaExplicit.h
 * @addtogroup Integrator
 * @{
 */

#ifndef GENERALIZEDALPHAEXPLICIT_H
#define GENERALIZEDALPHAEXPLICIT_H

#include "../Integrator.h"

class GeneralizedAlphaExplicit final : public ExplicitIntegrator {
    const double B, AM, AF;
    double DT{0.};

protected:
    void update_parameter(double) override;

    [[nodiscard]] int process_load_impl(bool) override;
    [[nodiscard]] int process_constraint_impl(bool) override;

    [[nodiscard]] bool has_corrector() const override;

    int correct_trial_status() override;

public:
    GeneralizedAlphaExplicit(unsigned, double);

    void assemble_resistance() override;

    vec get_force_residual() override;
    vec get_displacement_residual() override;
    sp_mat get_reference_load() override;

    int update_trial_status(bool) override;

    vec from_incre_acceleration(const vec&, const uvec&) override;
    vec from_total_acceleration(const vec&, const uvec&) override;

    void print() override;
};

#endif

//! @}
