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
 * @class NonlinearPeric
 * @brief The NonlinearPeric class.
 *
 * algorithm verified at 24 April 2019 by tlc
 *
 * @author tlc
 * @date 24/04/2019
 * @version 1.0.0
 * @file NonlinearPeric.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NONLINEARPERIC_H
#define NONLINEARPERIC_H

#include <Material/Material3D/Material3D.h>

struct DataNonlinearPeric {
    const double elastic_modulus; // elastic modulus
    const double poissons_ratio;  // poisson's ratio
    const double mu;
    const double epsilon;
};

class NonlinearPeric : DataNonlinearPeric, public Material3D {
    static constexpr unsigned max_iteration = 20;
    static const double root_three_two;
    static const mat unit_dev_tensor;
    static const double unit_time;

    const double shear_modulus = elastic_modulus / (2. + 2. * poissons_ratio); // shear modulus
    const double double_shear = 2. * shear_modulus;                            // double shear modulus
    const double triple_shear = 3. * shear_modulus;                            // double shear modulus

    const double factor_a = mu * epsilon;

    const double* incre_time = nullptr;

    [[nodiscard]] virtual double compute_k(double) const = 0;
    [[nodiscard]] virtual double compute_dk(double) const = 0;

public:
    NonlinearPeric(unsigned,    // tag
                   double,      // elastic modulus
                   double,      // poisson's ratio
                   double = 0., // mu
                   double = 0., // epsilon
                   double = 0.  // density
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
