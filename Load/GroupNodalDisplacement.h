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
 * @class GroupNodalDisplacement
 * @brief A GroupNodalDisplacement class.
 *
 * The GroupNodalDisplacement class is in charge of handling displacement load.
 *
 * @author tlc
 * @date 20/03/2020
 * @version 0.1.0
 * @file GroupNodalDisplacement.h
 * @addtogroup Load
 * @{
 */

#ifndef GROUPNODALDISPLACEMENT_H
#define GROUPNODALDISPLACEMENT_H

#include "NodalDisplacement.h"

class GroupNodalDisplacement final : protected GroupLoad, public NodalDisplacement {
public:
    GroupNodalDisplacement(
        unsigned, // tag
        double,   // magnitude
        uvec&&,   // group tags
        uvec&&,   // dof tags
        unsigned  // amplitude tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

#endif // GROUPNODALDISPLACEMENT_H

//! @}
