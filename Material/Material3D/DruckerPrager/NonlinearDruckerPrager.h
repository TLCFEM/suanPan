/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class NonlinearDruckerPrager
 * @brief The NonlinearDruckerPrager class.
 *
 * algorithm verified at 24 April 2019 by tlc
 *
 * @author tlc
 * @date 24/04/2019
 * @version 1.0.0
 * @file NonlinearDruckerPrager.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NONLINEARDRUCKERPRAGER_H
#define NONLINEARDRUCKERPRAGER_H

#include <Material/Material3D/Material3D.h>
#include <Toolbox/utility.h>

struct DataNonlinearDruckerPrager {
    const double elastic_modulus; // elastic modulus
    const double poissons_ratio;  // poisson's ratio
    const double eta_yield;
    const double eta_flow;
    const double xi;
};

class NonlinearDruckerPrager : protected DataNonlinearDruckerPrager, public Material3D {
    static constexpr unsigned max_iteration = 20u;
    static const mat unit_dev_tensor;
    static const mat unit_x_unit;

    const double shear = elastic_modulus / (2. + 2. * poissons_ratio); // shear modulus
    const double bulk = elastic_modulus / (3. - 6. * poissons_ratio);  // bulk modulus
    const double double_shear = 2. * shear;                            // double shear modulus

    const double factor_a = shear + bulk * eta_flow * eta_yield;
    const double factor_b = xi * xi / eta_flow / eta_yield;
    const double factor_c = sqrt(2.) * shear * bulk;
    const double factor_d = bulk * bulk * eta_flow * eta_yield;

    const bool associated = suanpan::approx_equal(eta_yield, eta_flow);

    [[nodiscard]] virtual double compute_c(double) const = 0;
    [[nodiscard]] virtual double compute_dc(double) const = 0;

public:
    NonlinearDruckerPrager(
        unsigned,   // tag
        double,     // elastic modulus
        double,     // poisson's ratio
        double,     // eta_yield (hydrostatic stress related)
        double,     // eta_flow (dilatancy angle related)
        double,     // xi (cohesion related)
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
