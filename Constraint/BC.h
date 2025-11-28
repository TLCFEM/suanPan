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
 * @class PenaltyBC
 * @class MultiplierBC
 * @class GroupPenaltyBC
 * @class GroupMultiplierBC
 * @brief A PenaltyBC class handles boundary conditions.
 *
 * The PenaltyBC class is in charge of applying boundary conditions to the system. The
 * PenaltyBC class only takes care of homogeneous Dirichlet conditions. Non-homogeneous
 * displacement boundary conditions are treated as Load so that can be solved
 * iteratively. Others are handled by general constraint class such as MPC. The
 * PenaltyBC class stores the boundary condition category, type, node(s) and
 * corresponding DoF(s). The Domain invokes `process(const shared_ptr<Domain>&)`
 * method to modify the global stiffness matrix.
 *
 * @author tlc
 * @date 26/11/2025
 * @version 0.2.0
 * @file BC.h
 * @addtogroup Constraint
 * @{
 */

#ifndef BC_H
#define BC_H

#include "Constraint.h"

class PenaltyBC : public Constraint {
    static double multiplier;

    friend void set_constraint_multiplier(double);

    [[nodiscard]] bool validate_node() const final { return true; }

public:
    PenaltyBC(
        unsigned,
        uvec&&,                  // node tags
        std::vector<Node::DOF>&& // dof components
    );

    int process(const shared_ptr<DomainBase>&) override;
    int process_resistance(const shared_ptr<DomainBase>&) final;
};

class MultiplierBC : public PenaltyBC {
public:
    using PenaltyBC::PenaltyBC;

    int process(const shared_ptr<DomainBase>&) override;
};

class GroupPenaltyBC : protected GroupModifier, public MultiplierBC {
public:
    GroupPenaltyBC(
        unsigned,
        uvec&&,                  // node group tags
        std::vector<Node::DOF>&& // dof components
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) override;
};

class GroupMultiplierBC final : public GroupPenaltyBC {
public:
    using GroupPenaltyBC::GroupPenaltyBC;

    int process(const shared_ptr<DomainBase>&) override;
};

void set_constraint_multiplier(double);

#endif

//! @}
