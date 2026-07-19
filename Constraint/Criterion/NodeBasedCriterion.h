/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
 * @class NodeBasedCriterion
 * @brief A NodeBasedCriterion class.
 *
 * The NodeBasedCriterion class.
 *
 * @author tlc
 * @date 18/09/2020
 * @version 0.1.0
 * @file NodeBasedCriterion.h
 * @addtogroup Criterion
 * @{
 */

#ifndef NODEBASEDCRITERION_H
#define NODEBASEDCRITERION_H

#include <Constraint/Criterion/Criterion.h>
#include <Domain/Node.h>

class NodeBasedCriterion : public Criterion {
protected:
    const double limit;
    const unsigned node;
    const Node::DOF dof;

    std::vector<uword> get_dof(const shared_ptr<DomainBase>&) const;

public:
    NodeBasedCriterion(
        unsigned,  // tag
        unsigned,  // step tag
        unsigned,  // node tag
        Node::DOF, // dof tag
        double     // limit
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
