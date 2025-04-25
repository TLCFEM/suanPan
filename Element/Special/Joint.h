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
 * @class Joint
 * @brief The Joint class.
 *
 * @author tlc
 * @date 20/10/2020
 * @version 0.1.0
 * @file Joint.h
 * @addtogroup Special
 * @ingroup Element
 * @{
 */

#ifndef JOINT_H
#define JOINT_H

#include <Element/MaterialElement.h>

class Joint final : public MaterialElement1D {
    static constexpr unsigned j_node = 2;

    const unsigned j_dof;

    const unsigned j_size = j_node * j_dof;

    std::vector<unique_ptr<Material>> j_material;

public:
    Joint(
        unsigned, // tag
        uvec&&,   // node tags
        uvec&&    // material tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    std::vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
