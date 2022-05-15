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
 * @class Damper02
 * @brief A Damper02 class.
 * @author tlc
 * @date 19/10/2017
 * @file Damper02.h
 * @addtogroup Special
 * @ingroup Element
 * @{
 */

#ifndef DAMPER02_H
#define DAMPER02_H

#include <Element/MaterialElement.h>

class Damper02 final : public MaterialElement1D {
    static constexpr unsigned d_node = 2, d_dof = 2, d_size = d_dof * d_node;

    static uvec IS, JS;

    const vec direction_cosine;

    unique_ptr<Material> device;

public:
    Damper02(unsigned, // tag
             uvec&&,   // node tag
             unsigned, // damper tag
             unsigned, // spring tag
             bool,     // if to use matrix formulation
             unsigned, // if proceed when fail to converge
             double    // beta
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
