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
 * @class CustomElastic1D
 * @brief A 1D Elastic class using custom constitutive equation.
 * @author tlc
 * @date 12/01/2023
 * @file CustomElastic1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CUSTOMELASTIC1D_H
#define CUSTOMELASTIC1D_H

#include <Material/Material1D/Material1D.h>
#include <Toolbox/Expression.h>
#include <Toolbox/ResourceHolder.h>

class CustomElastic1D final : public Material1D {
    const unsigned expression_tag;

    ResourceHolder<Expression> expression;

public:
    CustomElastic1D(
        unsigned,   // tag
        unsigned,   // expression tag
        double = 0. // density
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
