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
 * @class StressWrapper
 * @brief A StressWrapper class.
 *
 *  A universal wrapper that wraps 3D material models into other models by setting some stress components to zero.
 *
 * @author tlc
 * @date 22/09/2023
 * @version 0.1.0
 * @file StressWrapper.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef StressWrapper_H
#define StressWrapper_H

#include <Material/Material.h>
#include <Toolbox/ResourceHolder.h>

class StressWrapper : public Material {
    const uvec F1, F2;

    const unsigned base_tag;

    const unsigned max_iteration;

    vec trial_full_strain, current_full_strain;

    [[nodiscard]] mat form_stiffness(const mat&) const;

protected:
    ResourceHolder<Material> base;

public:
    StressWrapper(
        unsigned,    // tag
        unsigned,    // 3D material tag
        unsigned,    // max iteration
        uvec&&,      // non-trivial stress DoF
        uvec&&,      // trivial stress DoF
        MaterialType // material type
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    std::vector<vec> record(OutputType) override;
};

#endif

//! @}
