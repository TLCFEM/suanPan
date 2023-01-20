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
 * @class GeneralizedAlpha
 * @brief A GeneralizedAlpha class defines a solver using GeneralizedAlpha
 * algorithm.
 *
 * Unlike Newmark method, in which the equilibrium is satisfied at the end of
 * current time step, i.e., \f$t=t_0+\Delta{}t\f$, the generalized-\f$\alpha\f$
 * approach applies it at somewhere in current step, i.e.,
 * \f$t=t_0+\Delta{}t-\alpha\f$, similar to the generalized midpoint concept.
 *
 * doi: 10.1115/1.2900803
 *
 * @author tlc
 * @date 21/10/2017
 * @version 0.1.0
 * @file GeneralizedAlpha.h
 * @addtogroup Integrator
 * @{
 */

#ifndef GENERALIZEDALPHA_H
#define GENERALIZEDALPHA_H

#include "Integrator.h"

class GeneralizedAlpha final : public ImplicitIntegrator {
    const double alpha_f;
    const double alpha_m;
    const double gamma;
    const double beta;

    const double F1, F2, F3, F4, F9;

    double F5 = 0., F6 = 0., F7 = 0., F8 = 0., F10 = 0., F11 = 0.;

public:
    GeneralizedAlpha(unsigned, double);
    GeneralizedAlpha(unsigned, double, double);

    void assemble_resistance() override;
    void assemble_matrix() override;

    vec get_force_residual() override;
    vec get_displacement_residual() override;
    sp_mat get_reference_load() override;

    [[nodiscard]] int process_load() override;
    [[nodiscard]] int process_constraint() override;
    [[nodiscard]] int process_load_resistance() override;
    [[nodiscard]] int process_constraint_resistance() override;

    int update_trial_status() override;

    void update_parameter(double) override;

    vec from_incre_velocity(const vec&, const uvec&) override;
    vec from_incre_acceleration(const vec&, const uvec&) override;

    void print() override;
};

#endif

//! @}
