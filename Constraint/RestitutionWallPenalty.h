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
#include <Toolbox/container.h>

class Node;

class RestitutionWallPenalty : public RigidWallPenalty {
    suanpan::set<shared_ptr<Node>> node_pool;

    const double restitution_coefficient;

public:
    RestitutionWallPenalty(unsigned, unsigned, unsigned, vec&&, vec&&, double, double, unsigned);
    RestitutionWallPenalty(unsigned, unsigned, unsigned, vec&&, vec&&, vec&&, double, double, unsigned);

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) override;

    void stage(const shared_ptr<DomainBase>&) override;

    void commit_status() override;
    void clear_status() override;
    void reset_status() override;
};

class RestitutionWallPenalty1D final : public RestitutionWallPenalty {
public:
    RestitutionWallPenalty1D(unsigned, unsigned, unsigned, vec&&, vec&&, double, double);
};

class RestitutionWallPenalty2D final : public RestitutionWallPenalty {
public:
    RestitutionWallPenalty2D(unsigned, unsigned, unsigned, vec&&, vec&&, double, double);
    RestitutionWallPenalty2D(unsigned, unsigned, unsigned, vec&&, vec&&, vec&&, double, double);
};

class RestitutionWallPenalty3D final : public RestitutionWallPenalty {
public:
    RestitutionWallPenalty3D(unsigned, unsigned, unsigned, vec&&, vec&&, double, double);
    RestitutionWallPenalty3D(unsigned, unsigned, unsigned, vec&&, vec&&, vec&&, double, double);
};

#endif

//! @}
