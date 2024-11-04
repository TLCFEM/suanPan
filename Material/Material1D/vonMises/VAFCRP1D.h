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
 * @class VAFCRP1D
 * @brief The VAFCRP1D class defines a nonlinear hardening material with mixed
 * hardening (isotropic and kinematic) based on J2 plasticity rule.
 *
 * The isotropic hardening is defined as an exponential function.
 *
 * The kinematic hardening consists of multiple Armstrong--Frederick type back stresses.
 *
 * algorithm verified at 28 October 2019 by tlc
 *
 * @author tlc
 * @date 03/11/2024
 * @version 2.0.0
 * @file VAFCRP1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef VAFCRP1D_H
#define VAFCRP1D_H

#include <Material/Material1D/Material1D.h>

struct DataVAFCRP1D {
    const double elastic_modulus; // elastic modulus
    const double yield;           // yield stress
    const double hardening;       // linear isotropic hardening modulus
    const double saturated;       // saturation stress
    const double m;               // saturation rate
    const double mu, epsilon;
    const vec a, b;
};

class VAFCRP1D final : protected DataVAFCRP1D, public Material1D {
    static constexpr unsigned max_iteration = 20u;

    const double* incre_time = nullptr;

    const unsigned size = static_cast<unsigned>(a.size());

    [[nodiscard]] std::tuple<double, double> isotropic_bound(double) const;

    [[nodiscard]] bool is_elastic(double, const vec&) const;

    int partial_loading(const vec&, double);

public:
    VAFCRP1D(
        unsigned,       // tag
        DataVAFCRP1D&&, // material data
        double = 0.     // density
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
