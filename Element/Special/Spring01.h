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
 * @class Spring01
 * @brief The Spring01 class.
 *
 * @author tlc
 * @date 02/03/2019
 * @version 0.2.0
 * @file Spring01.h
 * @addtogroup Special
 * @ingroup Element
 * @{
 */

#ifndef SPRING01_H
#define SPRING01_H

#include <Element/MaterialElement.h>

class Spring01 final : public MaterialElement1D {
    static constexpr unsigned s_node = 2, s_dof = 2, s_size = s_dof * s_node;

    static uvec IS, JS;

    vec direction_cosine;

    unique_ptr<Material> s_material;

public:
    Spring01(
        unsigned, // tag
        uvec&&,   // node tags
        unsigned  // material tag
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
