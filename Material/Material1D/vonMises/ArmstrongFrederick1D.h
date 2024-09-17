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
 * @class ArmstrongFrederick1D
 * @brief The ArmstrongFrederick1D class defines a nonlinear hardening material with mixed
 * hardening (isotropic and kinematic) based on J2 plasticity rule.
 *
 * The isotropic hardening is defined as an exponential function.
 *
 * The kinematic hardening consists of multiple Armstrong--Frederick type back stresses.
 *
 * algorithm verified at 06 October 2019 by tlc
 *
 * @author tlc
 * @date 06/10/2019
 * @version 1.0.0
 * @file ArmstrongFrederick1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef ARMSTRONGFREDERICK1D_H
#define ARMSTRONGFREDERICK1D_H

#include <Material/Material1D/Material1D.h>

struct DataArmstrongFrederick1D {
    const double elastic_modulus; // elastic modulus
    const double yield;           // yield stress
    const double hardening;       // linear isotropic hardening modulus
    const double saturation;      // saturation stress
    const double ms;              // saturation rate
    const double memory;          // strain memory ratio
    const double reduction;       // isotropic hardening reduction
    const double mr;              // isotropic hardening reduction rate
    const vec a, b;
};

class ArmstrongFrederick1D final : protected DataArmstrongFrederick1D, public Material1D {
    static constexpr unsigned max_iteration = 20u;

    const unsigned size = static_cast<unsigned>(a.size());

public:
    ArmstrongFrederick1D(
        unsigned,                   // tag
        DataArmstrongFrederick1D&&, // material data
        double = 0.                 // density
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
