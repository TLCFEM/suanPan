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
 * @class GroupMultiplierBC
 * @brief A GroupMultiplierBC class handles boundary conditions.
 *
 * The GroupMultiplierBC class is in charge of applying boundary conditions to the system. The
 * GroupMultiplierBC class only takes care of homogeneous Dirichlet conditions. Non-homogeneous
 * displacement boundary conditions are treated as Load so that can be solved
 * iteratively. Others are handled by general constraint class such as MPC. The
 * GroupMultiplierBC class stores the boundary condition category, type, node(s) and
 * corresponding DoF(s). The Domain invokes `process(const shared_ptr<Domain>&)`
 * method to modify the global stiffness matrix.
 *
 * @author tlc
 * @date 21/03/2020
 * @version 0.1.0
 * @file GroupMultiplierBC.h
 * @addtogroup Constraint
 * @{
 */

#ifndef GROUPMULTIPLIERBC_H
#define GROUPMULTIPLIERBC_H

#include "GroupPenaltyBC.h"

class GroupMultiplierBC final : public GroupPenaltyBC {
public:
    using GroupPenaltyBC::GroupPenaltyBC;

    int process(const shared_ptr<DomainBase>&) override;
    int process_resistance(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
