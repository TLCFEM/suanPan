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
 * @class NonlinearViscosity
 * @brief A 1D Viscosity class.
 * @author tlc
 * @date 01/09/2020
 * @version 0.2.0
 * @file NonlinearViscosity.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef NONLINEARVISCOSITY_H
#define NONLINEARVISCOSITY_H

#include <Material/Material1D/Material1D.h>

struct DataNonlinearViscosity {
    const double alpha;
    const double limit; // cubic replacement
};

class NonlinearViscosity : protected DataNonlinearViscosity, public Material1D {
    const double a = (.5 * alpha - .5) * pow(limit, alpha - 3.);
    const double b = (1.5 - .5 * alpha) * pow(limit, alpha - 1.);

    [[nodiscard]] virtual double compute_du(double, double) const = 0;
    [[nodiscard]] virtual double compute_dv(double, double) const = 0;
    [[nodiscard]] virtual double compute_damping_coefficient(double, double) const = 0;

public:
    NonlinearViscosity(
        unsigned, // tag
        double,   // alpha
        double    // cut-off
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) final;
    int update_trial_status(const vec&, const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    std::vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
