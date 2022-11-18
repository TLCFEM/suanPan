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
 * @class BilinearJ2
 * @brief The BilinearJ2 class defines a bilinear hardening material with mixed
 * hardening (isotropic and kinematic) based on J2 plasticity rule.
 *
 * algorithm verified at 24 April 2019 by tlc
 *
 * @author tlc
 * @date 24/04/2019
 * @version 1.0.0
 * @file BilinearJ2.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef BILINEARJ2_H
#define BILINEARJ2_H

#include <Material/Material3D/Material3D.h>

struct DataBilinearJ2 {
    const double elastic_modulus; // elastic modulus
    const double poissons_ratio;  // poisson's ratio
    const double yield_stress;    // initial yield stress
    const double hardening_ratio; // hardening ratio
    const double beta;            // isotropic (1.0) / kinematic (0.0) hardening factor
};

class BilinearJ2 final : DataBilinearJ2, public Material3D {
    static const double two_third;
    static const double root_two_third;
    static const mat66 unit_dev_tensor;

    const double shear_modulus = elastic_modulus / (2. + 2. * poissons_ratio); // shear modulus
    const double double_shear = 2. * shear_modulus;                            // double shear modulus
    const double square_double_shear = double_shear * double_shear;            // double shear modulus squared
    const double isotropic_modulus = beta * elastic_modulus * hardening_ratio / (1. - hardening_ratio);
    const double kinematic_modulus = (1. - beta) * elastic_modulus * hardening_ratio / (1. - hardening_ratio);

public:
    BilinearJ2(unsigned,    // tag
               double,      // elastic modulus
               double,      // poisson's ratio
               double,      // initial yield stress
               double = 0., // hardening ratio
               double = 1., // isotropic (1.0) / kinematic (0.0) hardening factor
               double = 0.  // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

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
