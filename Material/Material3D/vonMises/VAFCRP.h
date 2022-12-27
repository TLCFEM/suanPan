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
 * @class VAFCRP
 * @brief The VAFCRP class defines a nonlinear hardening material with mixed
 * hardening (isotropic and kinematic) based on J2 plasticity rule.
 *
 * The isotropic hardening is defined as a Voce type exponential function.
 *
 * The kinematic hardening consists of multiple Armstrong--Frederick type back stresses.
 *
 * The Perzyna type viscosity is used.
 *
 * algorithm verified at 25 October 2019 by tlc
 *
 * @author tlc
 * @date 25/10/2019
 * @version 1.0.0
 * @file VAFCRP.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef VAFCRP_H
#define VAFCRP_H

#include <Material/Material3D/Material3D.h>

struct DataVAFCRP {
    const double elastic_modulus; // elastic modulus
    const double poissons_ratio;  // poisson's ratio
    const double yield;           // yield stress
    const double saturated;
    const double hardening;
    const double m;
    const double mu;
    const double epsilon;
    const vec a, b;
};

class VAFCRP final : DataVAFCRP, public Material3D {
    static constexpr unsigned max_iteration = 20;
    static const double root_three_two;
    static const mat unit_dev_tensor;

    const double* incre_time = nullptr;

    const unsigned size = static_cast<unsigned>(a.size());

    const double shear = elastic_modulus / (2. + 2. * poissons_ratio); // shear modulus
    const double double_shear = 2. * shear;
    const double three_shear = 3. * shear;
    const double root_six_shear = sqrt(6.) * shear;

public:
    VAFCRP(unsigned,   // tag
           double,     // elastic modulus
           double,     // poissons ratio
           double,     // yield stress
           double,     // saturated stress
           double,     // linear hardening modulus
           double,     // m
           double,     // mu
           double,     // epsilon
           vec&&,      // a
           vec&&,      // b
           double = 0. // density
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
