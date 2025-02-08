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
 * @class NodeGroup
 * @brief The NodeGroup class.
 *
 * @author tlc
 * @date 21/05/2020
 * @version 0.1.0
 * @file NodeGroup.h
 * @addtogroup Group
 * @{
 */

#ifndef NODEGROUP_H
#define NODEGROUP_H

#include "Group.h"

class NodeGroup final : public Group {
    const int dof;

    const vec rule;

    const vec s_node, e_node;

public:
    NodeGroup(unsigned, int, vec&&);
    NodeGroup(unsigned, uvec&&);
    NodeGroup(unsigned, vec&&, vec&&);
    NodeGroup(unsigned, vec&&);

    void initialize(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
