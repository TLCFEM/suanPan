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
 * @class SimpleSand
 * @brief The SimpleSand class.
 *
 * verified and approved by tlc @ 12/07/2021
 *
 * @author tlc
 * @date 12/07/2021
 * @version 0.1.0
 * @file SimpleSand.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef SIMPLESAND_H
#define SIMPLESAND_H

#include <Material/Material3D/Material3D.h>

struct DataSimpleSand {
    const double elastic_modulus; // elastic modulus
    const double poissons_ratio;  // poisson's ratio
    const double m = .01;
    const double a = -.7;
    const double h = 5.;
    const double ac = 1.25;
    const double nb = 1.1;
    const double nd = 3.5;
    const double vc = 1.915;
    const double pc = -130.;
    const double lc = .02;
    const double v0 = 2.;
};

class SimpleSand final : DataSimpleSand, public Material3D {
    static constexpr unsigned max_iteration = 20;
    static const mat unit_dev_tensor;

    static constexpr uword sa = 0, sb = 1;
    static const span sc, sd;

    const double shear = elastic_modulus / (2. + 2. * poissons_ratio); // shear modulus
    const double double_shear = 2. * shear;                            // double shear
    const double bulk = elastic_modulus / (3. - 6. * poissons_ratio);

public:
    SimpleSand(unsigned,   // tag
               double,     // elastic modulus
               double,     // poissons ratio
               double,     // m
               double,     // A
               double,     // h
               double,     // alpha_c
               double,     // n_b
               double,     // n_d
               double,     // v_c
               double,     // p_c
               double,     // lambda_c
               double,     // v_0
               double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
