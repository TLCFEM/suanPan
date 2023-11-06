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
 * @class Sequential
 * @brief A Sequential material class.
 *
 * The Sequential material class is a container to hold multiple 1D material
 * models together one next to another. For each material model, the strain may
 * be different but the stress is identical.
 *
 * @author tlc
 * @date 08/06/2018
 * @version 0.1.0
 * @file Sequential.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef SEQUENTIAL_H
#define SEQUENTIAL_H

#include <Material/Material1D/Material1D.h>
#include <Toolbox/ResourceHolder.h>

class Sequential final : public Material1D {
    static constexpr unsigned max_iteration = 20u;

    const uword mat_size;
    const uvec mat_tag;
    vector<ResourceHolder<Material>> mat_pool;
    mat jacobian;

public:
    Sequential(
        unsigned, // tag
        uvec&&    // material tag pool
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;
    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif
