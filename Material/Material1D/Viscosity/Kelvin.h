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
 * @class Kelvin
 * @brief A 1D Kelvin material class.
 * @author tlc
 * @date 13/01/2021
 * @version 0.2.0
 * @file KELVIN.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef KELVIN_H
#define KELVIN_H

#include <Material/Material1D/Material1D.h>
#include <Toolbox/ResourceHolder.h>

class Kelvin final : public Material1D {
    const unsigned damper_tag, spring_tag;

    ResourceHolder<Material> damper, spring;

public:
    Kelvin(
        unsigned, // tag
        unsigned, // damper tag
        unsigned  // spring tag
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

//! @}
