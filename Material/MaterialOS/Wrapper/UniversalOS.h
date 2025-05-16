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
 * @class UniversalOS
 * @brief A UniversalOS class.
 * @author tlc
 * @date 21/09/2023
 * @version 0.1.0
 * @file UniversalOS.h
 * @addtogroup Material-OS
 * @{
 */

#ifndef UNIVERSALOS_H
#define UNIVERSALOS_H

#include <Material/Material3D/Wrapper/StressWrapper.h>

class UniversalOS : public StressWrapper {
public:
    UniversalOS(
        unsigned, // tag
        unsigned, // 3D material tag
        unsigned, // max iteration
        uvec&&,
        uvec&&
    );

    void print() override;
};

class OS146 final : public UniversalOS {
public:
    OS146(
        unsigned, // tag
        unsigned, // 3D material tag
        unsigned  // max iteration
    );

    unique_ptr<Material> get_copy() override;
};

class OS146S final : public Material {
    const unsigned base_tag;

    const double shear_modulus;

    ResourceHolder<Material> base;

public:
    OS146S(
        unsigned, // tag
        unsigned, // 3D material tag
        double    // shear modulus
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    std::vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
