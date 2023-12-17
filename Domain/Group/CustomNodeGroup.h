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
 * @class CustomNodeGroup
 * @brief The CustomNodeGroup class.
 *
 * @author tlc
 * @date 04/03/2023
 * @version 0.1.0
 * @file CustomNodeGroup.h
 * @addtogroup Group
 * @{
 */

#ifndef CUSTOMNODEGROUP_H
#define CUSTOMNODEGROUP_H

#include "Group.h"
#include <Toolbox/Expression.h>
#include <Toolbox/ResourceHolder.h>

class CustomNodeGroup final : public Group {
    const unsigned expression_tag;

    ResourceHolder<Expression> expression;

public:
    CustomNodeGroup(unsigned, unsigned);

    void initialize(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
