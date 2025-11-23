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

#ifndef NODE_HELPER_HPP
#define NODE_HELPER_HPP

#include "Node.h"

template<Node::DOF...> struct always_false {
    static constexpr bool value = false;
};

template<Node::DOF... D> bool check_dof_definition(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_current_position(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_trial_position(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_current_displacement(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_current_velocity(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_current_acceleration(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_trial_displacement(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_trial_velocity(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_trial_acceleration(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_incre_displacement(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_incre_velocity(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<Node::DOF... D> vec get_incre_acceleration(const shared_ptr<Node>&) {
    static_assert(always_false<D...>::value, "not implemented");
    return {};
}

template<> inline bool check_dof_definition<Node::DOF::U1>(const shared_ptr<Node>& t_node) {
    auto& t_dof = t_node->get_dof_identifier();
    return !t_dof.empty() && t_dof.at(0) == Node::DOF::U1;
}

template<> inline bool check_dof_definition<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) {
    auto& t_dof = t_node->get_dof_identifier();
    return t_dof.size() > 1 && t_dof.at(0) == Node::DOF::U1 && t_dof.at(1) == Node::DOF::U2;
}

template<> inline bool check_dof_definition<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) {
    auto& t_dof = t_node->get_dof_identifier();
    return t_dof.size() > 2 && t_dof.at(0) == Node::DOF::U1 && t_dof.at(1) == Node::DOF::U2 && t_dof.at(2) == Node::DOF::U3;
}

template<> inline vec get_current_position<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_displacement()).resize(1) + vec(t_node->get_coordinate()).resize(1); }

template<> inline vec get_trial_position<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_displacement()).resize(1) + vec(t_node->get_coordinate()).resize(1); }

template<> inline vec get_current_displacement<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_displacement()).resize(1); }

template<> inline vec get_current_velocity<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_velocity()).resize(1); }

template<> inline vec get_current_acceleration<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_acceleration()).resize(1); }

template<> inline vec get_trial_displacement<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_displacement()).resize(1); }

template<> inline vec get_trial_velocity<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_velocity()).resize(1); }

template<> inline vec get_trial_acceleration<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_acceleration()).resize(1); }

template<> inline vec get_incre_displacement<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_incre_displacement()).resize(1); }

template<> inline vec get_incre_velocity<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_incre_velocity()).resize(1); }

template<> inline vec get_incre_acceleration<Node::DOF::U1>(const shared_ptr<Node>& t_node) { return vec(t_node->get_incre_acceleration()).resize(1); }

template<> inline vec get_current_position<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_displacement()).resize(2) + vec(t_node->get_coordinate()).resize(2); }

template<> inline vec get_trial_position<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_displacement()).resize(2) + vec(t_node->get_coordinate()).resize(2); }

template<> inline vec get_current_displacement<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_displacement()).resize(2); }

template<> inline vec get_current_velocity<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_velocity()).resize(2); }

template<> inline vec get_current_acceleration<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_acceleration()).resize(2); }

template<> inline vec get_trial_displacement<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_displacement()).resize(2); }

template<> inline vec get_trial_velocity<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_velocity()).resize(2); }

template<> inline vec get_trial_acceleration<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_acceleration()).resize(2); }

template<> inline vec get_incre_displacement<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_incre_displacement()).resize(2); }

template<> inline vec get_incre_velocity<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_incre_velocity()).resize(2); }

template<> inline vec get_incre_acceleration<Node::DOF::U1, Node::DOF::U2>(const shared_ptr<Node>& t_node) { return vec(t_node->get_incre_acceleration()).resize(2); }

template<> inline vec get_current_position<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_displacement()).resize(3) + vec(t_node->get_coordinate()).resize(3); }

template<> inline vec get_trial_position<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_displacement()).resize(3) + vec(t_node->get_coordinate()).resize(3); }

template<> inline vec get_current_displacement<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_displacement()).resize(3); }

template<> inline vec get_current_velocity<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_velocity()).resize(3); }

template<> inline vec get_current_acceleration<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_current_acceleration()).resize(3); }

template<> inline vec get_trial_displacement<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_displacement()).resize(3); }

template<> inline vec get_trial_velocity<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_velocity()).resize(3); }

template<> inline vec get_trial_acceleration<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_trial_acceleration()).resize(3); }

template<> inline vec get_incre_displacement<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_incre_displacement()).resize(3); }

template<> inline vec get_incre_velocity<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_incre_velocity()).resize(3); }

template<> inline vec get_incre_acceleration<Node::DOF::U1, Node::DOF::U2, Node::DOF::U3>(const shared_ptr<Node>& t_node) { return vec(t_node->get_incre_acceleration()).resize(3); }

#endif
