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
 * @class Damper01
 * @brief A Damper01 class.
 *
 * Using quadrant damper and displacement and velocity as basic input.
 *
 * @author tlc
 * @date 23/07/2018
 * @file Damper01.h
 * @addtogroup Special
 * @ingroup Element
 * @{
 */

#ifndef DAMPER01_H
#define DAMPER01_H

#include <Element/MaterialElement.h>

class Damper01 : public MaterialElement1D {
    static constexpr unsigned d_node = 2;

    const unsigned d_dof;

protected:
    const unsigned d_size = d_dof * d_node;

    const uvec IS, JS;

    const vec direction_cosine;

    unique_ptr<Material> damper;

public:
    Damper01(
        unsigned, // tag
        uvec&&,   // node tag
        unsigned, // damper tag
        unsigned  // dimension
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

class Damper05 final : public Damper01 {
public:
    using Damper01::Damper01;

    int update_status() override;
};

#endif

//! @}
