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
 * @class RigidWallPenalty
 * @brief A RigidWall class.
 *
 * @author tlc
 * @date 15/07/2020
 * @version 0.1.0
 * @file RigidWallPenalty.h
 * @addtogroup Constraint
 * @{
 */

#ifndef RIGIDWALLPENALTY_H
#define RIGIDWALLPENALTY_H

#include "Constraint.h"

class RigidWallPenalty : public Constraint {
    static auto setup(const unsigned dimension) {
        std::vector<Node::DOF> ref_dof;
        if(dimension > 0u) ref_dof.emplace_back(Node::DOF::U1);
        if(dimension > 1u) ref_dof.emplace_back(Node::DOF::U2);
        if(dimension > 2u) ref_dof.emplace_back(Node::DOF::U3);
        return ref_dof;
    }

    [[nodiscard]] bool reject_invalid_object() const override { return true; }

protected:
    const unsigned n_dim;

    const double alpha;

    const vec edge_a, edge_b;
    const vec origin, outer_norm;
    const double length_a = 0., length_b = 0.;

public:
    RigidWallPenalty(unsigned, unsigned, vec&&, vec&&, double, unsigned);
    RigidWallPenalty(unsigned, unsigned, vec&&, vec&&, vec&&, double, unsigned);

    int process(const shared_ptr<DomainBase>&) override;

    void commit_status() override;
    void clear_status() override;
    void reset_status() override;
};

class RigidWallPenalty1D final : public RigidWallPenalty {
public:
    RigidWallPenalty1D(unsigned, unsigned, double, double, double);
};

class RigidWallPenalty2D final : public RigidWallPenalty {
public:
    RigidWallPenalty2D(unsigned, unsigned, vec2&&, vec2&&, double);
    RigidWallPenalty2D(unsigned, unsigned, vec2&&, vec3&&, double);
};

class RigidWallPenalty3D final : public RigidWallPenalty {
public:
    RigidWallPenalty3D(unsigned, unsigned, vec3&&, vec3&&, double);
    RigidWallPenalty3D(unsigned, unsigned, vec3&&, vec3&&, vec3&&, double);
};

#endif

//! @}
