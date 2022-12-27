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
 * @class MultiplierBC
 * @brief A MultiplierBC class handles boundary conditions.
 *
 * The MultiplierBC class is in charge of applying boundary conditions to the system. The
 * MultiplierBC class only takes care of homogeneous Dirichlet conditions. Non-homogeneous
 * displacement boundary conditions are treated as Load so that can be solved
 * iteratively. Others are handled by general constraint class such as MPC. The
 * MultiplierBC class stores the boundary condition category, type, node(s) and
 * corresponding DoF(s). The Domain invokes `process(const shared_ptr<Domain>&)`
 * method to modify the global stiffness matrix.
 *
 * @author tlc
 * @date 23/07/2017
 * @version 0.1.0
 * @file MultiplierBC.h
 * @addtogroup Constraint
 * @{
 */

#ifndef MULTIPLIERBC_H
#define MULTIPLIERBC_H

#include "PenaltyBC.h"

class MultiplierBC : public PenaltyBC {
public:
    using PenaltyBC::PenaltyBC;

    int process(const shared_ptr<DomainBase>&) override;
    int process_resistance(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
