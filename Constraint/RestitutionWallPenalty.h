/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
 * @class RestitutionWallPenalty
 * @brief A RestitutionWallPenalty class.
 *
 * @author tlc
 * @date 04/06/2022
 * @version 0.1.0
 * @file RestitutionWallPenalty.h
 * @addtogroup Constraint
 * @{
 */

#ifndef RESTITUTIONWALLPENALTY_H
#define RESTITUTIONWALLPENALTY_H

#include "RigidWallPenalty.h"
#include "Toolbox/container.h"

class Node;

class RestitutionWallPenalty final : public RigidWallPenalty {
    suanpan::set<shared_ptr<Node>> node_pool;

    const double restitution_coefficient;

public:
    RestitutionWallPenalty(unsigned, unsigned, unsigned, vec&&, vec&&, double, double);
    RestitutionWallPenalty(unsigned, unsigned, unsigned, vec&&, vec&&, vec&&, double, double);

    int process(const shared_ptr<DomainBase>&) override;

    void commit_status() override;
    void clear_status() override;
    void reset_status() override;
};

#endif

//! @}
