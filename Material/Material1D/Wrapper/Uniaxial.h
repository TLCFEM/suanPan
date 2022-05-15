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
 * @class Uniaxial
 * @brief A Uniaxial class.
 * @author tlc
 * @date 22/01/2019
 * @version 0.1.0
 * @file Uniaxial.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef UNIAXIAL_H
#define UNIAXIAL_H

#include <Material/Material1D/Material1D.h>

class Uniaxial final : public Material1D {
    static const uvec F1, F2;

    const unsigned base_tag;

    const unsigned max_iteration;

    unique_ptr<Material> base;

    vec trial_full_strain, current_full_strain;

    static mat form_stiffness(const mat&);

public:
    Uniaxial(unsigned,    // tag
             unsigned,    // 3D material tag
             unsigned = 1 // max iteration
    );
    Uniaxial(const Uniaxial&);
    Uniaxial(Uniaxial&&) = delete;
    Uniaxial& operator=(const Uniaxial&) = delete;
    Uniaxial& operator=(Uniaxial&&) = delete;
    ~Uniaxial() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
