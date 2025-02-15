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
 * @class NonlinearGurson
 * @brief The NonlinearGurson class.
 *
 * @author tlc
 * @date 05/01/2020
 * @version 1.0.0
 * @file NonlinearGurson.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NONLINEARGURSON_H
#define NONLINEARGURSON_H

#include <Material/Material3D/Material3D.h>

struct DataNonlinearGurson {
    const double elastic_modulus; // elastic modulus
    const double poissons_ratio;  // poisson's ratio

    const double q1 = 1., q2 = 1., fn = 0., sn = 1., en = 0.;
};

class NonlinearGurson : protected DataNonlinearGurson, public Material3D {
    static constexpr unsigned max_iteration = 20u;
    static const double sqrt_three_two;
    static const mat unit_dev_tensor;

    const double six_shear = 3. * elastic_modulus / (1. + poissons_ratio); // six shear modulus
    const double bulk = elastic_modulus / (3. - 6. * poissons_ratio);      // bulk modulus

    const double para_a = 3. * bulk * q1 * q2;
    const double para_b = fn / sn / sqrt(2. * datum::pi);

    [[nodiscard]] virtual vec2 compute_hardening(double) const = 0;

public:
    NonlinearGurson(
        unsigned,   // tag
        double,     // elastic modulus
        double,     // poisson's ratio
        double,     // q1
        double,     // q2
        double,     // fn
        double,     // sn
        double,     // en
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
