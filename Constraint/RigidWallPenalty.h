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
#include <Domain/NodeHelper.hpp>

class RigidWallPenalty : public Constraint {
protected:
    template<DOF... D> void set_handler() { throw std::logic_error("not implemented"); }

    const unsigned n_dim;

    const double alpha;

    const vec edge_a, edge_b;
    const vec origin, outer_norm;
    const double length_a = 0., length_b = 0.;

    bool (*checker_handler)(const shared_ptr<Node>&) = nullptr;
    Col<double> (*current_velocity_handler)(const shared_ptr<Node>&) = nullptr;
    Col<double> (*incre_acceleration_handler)(const shared_ptr<Node>&) = nullptr;
    Col<double> (*trial_position_handler)(const shared_ptr<Node>&) = nullptr;
    Col<double> (*trial_displacement_handler)(const shared_ptr<Node>&) = nullptr;
    Col<double> (*trial_velocity_handler)(const shared_ptr<Node>&) = nullptr;
    Col<double> (*trial_acceleration_handler)(const shared_ptr<Node>&) = nullptr;

public:
    RigidWallPenalty(unsigned, unsigned, unsigned, vec&&, vec&&, double, unsigned);
    RigidWallPenalty(unsigned, unsigned, unsigned, vec&&, vec&&, vec&&, double, unsigned);

    int process(const shared_ptr<DomainBase>&) override;

    void commit_status() override;
    void clear_status() override;
    void reset_status() override;
};

template<> inline void RigidWallPenalty::set_handler<DOF::U1>() {
    checker_handler = check_dof_definition<DOF::U1>;
    current_velocity_handler = get_current_velocity<DOF::U1>;
    incre_acceleration_handler = get_incre_acceleration<DOF::U1>;
    trial_position_handler = get_trial_position<DOF::U1>;
    trial_displacement_handler = get_trial_displacement<DOF::U1>;
    trial_velocity_handler = get_trial_velocity<DOF::U1>;
    trial_acceleration_handler = get_trial_acceleration<DOF::U1>;
}

template<> inline void RigidWallPenalty::set_handler<DOF::U1, DOF::U2>() {
    checker_handler = check_dof_definition<DOF::U1, DOF::U2>;
    current_velocity_handler = get_current_velocity<DOF::U1, DOF::U2>;
    incre_acceleration_handler = get_incre_acceleration<DOF::U1, DOF::U2>;
    trial_position_handler = get_trial_position<DOF::U1, DOF::U2>;
    trial_displacement_handler = get_trial_displacement<DOF::U1, DOF::U2>;
    trial_velocity_handler = get_trial_velocity<DOF::U1, DOF::U2>;
    trial_acceleration_handler = get_trial_acceleration<DOF::U1, DOF::U2>;
}

template<> inline void RigidWallPenalty::set_handler<DOF::U1, DOF::U2, DOF::U3>() {
    checker_handler = check_dof_definition<DOF::U1, DOF::U2, DOF::U3>;
    current_velocity_handler = get_current_velocity<DOF::U1, DOF::U2, DOF::U3>;
    incre_acceleration_handler = get_incre_acceleration<DOF::U1, DOF::U2, DOF::U3>;
    trial_position_handler = get_trial_position<DOF::U1, DOF::U2, DOF::U3>;
    trial_displacement_handler = get_trial_displacement<DOF::U1, DOF::U2, DOF::U3>;
    trial_velocity_handler = get_trial_velocity<DOF::U1, DOF::U2, DOF::U3>;
    trial_acceleration_handler = get_trial_acceleration<DOF::U1, DOF::U2, DOF::U3>;
}

class RigidWallPenalty1D final : public RigidWallPenalty {
public:
    RigidWallPenalty1D(unsigned, unsigned, unsigned, vec&&, vec&&, double);
};

class RigidWallPenalty2D final : public RigidWallPenalty {
public:
    RigidWallPenalty2D(unsigned, unsigned, unsigned, vec&&, vec&&, double);
    RigidWallPenalty2D(unsigned, unsigned, unsigned, vec&&, vec&&, vec&&, double);
};

class RigidWallPenalty3D final : public RigidWallPenalty {
public:
    RigidWallPenalty3D(unsigned, unsigned, unsigned, vec&&, vec&&, double);
    RigidWallPenalty3D(unsigned, unsigned, unsigned, vec&&, vec&&, vec&&, double);
};

#endif

//! @}
