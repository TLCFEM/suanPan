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
 * @class Stacked
 * @brief A Stacked material class.
 * @author tlc
 * @date 20/02/2019
 * @version 0.1.1
 * @file Stacked.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef STACKED_H
#define STACKED_H

#include <Material/Material3D/Material3D.h>
#include <Toolbox/ResourceHolder.h>

class Stacked final : public Material3D {
    const uvec mat_tag;

    std::vector<ResourceHolder<Material>> mat_pool;

public:
    Stacked(
        unsigned, // tag
        uvec&&    // mat tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;

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
