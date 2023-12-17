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
 * @class SteelBRB
 * @brief The SteelBRB class.
 *
 * The BRB steel model.
 *
 * @author tlc
 * @date 19/08/2020
 * @version 1.0.0
 * @file SteelBRB.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef STEELBRB_H
#define STEELBRB_H

#include <Material/Material1D/Material1D.h>

struct DataSteelBRB {
    const double elastic_modulus = 2E5; // elastic modulus
    const double yield_stress = 400.;   // yield stress
    const double plastic_modulus = 2E3;
    const double t_saturated_stress = 660.;
    const double t_scalar = .2;
    const double t_exponent = .6;
    const double c_saturated_stress = 450.;
    const double c_scalar = .15;
    const double c_exponent = .4;
};

class SteelBRB final : protected DataSteelBRB, public Material1D {
    static constexpr unsigned max_iteration = 20u;

    const double s_modulus = -elastic_modulus - plastic_modulus;
    const double c_const = (c_saturated_stress - yield_stress) / c_scalar;
    const double t_const = (t_saturated_stress - yield_stress) / t_scalar;

    [[nodiscard]] vec compute_t_yield_stress(double) const;
    [[nodiscard]] vec compute_c_yield_stress(double) const;

public:
    SteelBRB(
        unsigned, // tag
        vec&&     // parameter
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
