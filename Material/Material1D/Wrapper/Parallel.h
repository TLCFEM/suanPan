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
 * @class Parallel
 * @brief A Parallel material class.
 *
 * The Parallel material class is a container to hold multiple 1D material models
 * together side by side. The strain is identical for all models and the overall
 * stress/stiffness is the summation of all stress/stiffness components.
 *
 * @author tlc
 * @date 08/06/2018
 * @version 0.1.0
 * @file Parallel.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef PARALLEL_H
#define PARALLEL_H

#include <Material/Material1D/Material1D.h>
#include <Toolbox/ResourceHolder.h>

class Parallel final : public Material1D {
    const uvec mat_tag;

    vector<ResourceHolder<Material>> mat_pool;

public:
    Parallel(
        unsigned, // tag
        uvec&&    // material tag pool
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;
    int update_trial_status(const vec&, const vec&) override;
    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif
