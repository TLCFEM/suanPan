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

#ifndef NODE_HELPER_HPP
#define NODE_HELPER_HPP

#include "DOF.h"
#include "Node.h"

template<DOF... D> bool check_dof_definition(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_current_position(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_trial_position(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_current_displacement(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_current_velocity(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_current_acceleration(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_trial_displacement(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_trial_velocity(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_trial_acceleration(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_incre_displacement(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_incre_velocity(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<DOF... D> vec get_incre_acceleration(const shared_ptr<Node>&) { throw std::logic_error("not implemented"); }

template<> inline bool check_dof_definition<DOF::U1>(const shared_ptr<Node>& t_node) {
    auto& t_dof = t_node->get_dof_identifier();
    return !t_dof.empty() && t_dof.at(0) == DOF::U1;
}

template<> inline bool check_dof_definition<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) {
    auto& t_dof = t_node->get_dof_identifier();
    return t_dof.size() > 1 && t_dof.at(0) == DOF::U1 && t_dof.at(1) == DOF::U2;
}

template<> inline bool check_dof_definition<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) {
    auto& t_dof = t_node->get_dof_identifier();
    return t_dof.size() > 2 && t_dof.at(0) == DOF::U1 && t_dof.at(1) == DOF::U2 && t_dof.at(2) == DOF::U3;
}

template<> inline vec get_current_position<DOF::U1>(const shared_ptr<Node>& t_node) {
    vec t_vec = t_node->get_current_displacement().head(1);
    if(auto& t_coor = t_node->get_coordinate(); !t_coor.empty()) t_vec(0) += t_coor(0);
    return t_vec;
}

template<> inline vec get_trial_position<DOF::U1>(const shared_ptr<Node>& t_node) {
    vec t_vec = t_node->get_trial_displacement().head(1);
    if(auto& t_coor = t_node->get_coordinate(); !t_coor.empty()) t_vec(0) += t_coor(0);
    return t_vec;
}

template<> inline vec get_current_displacement<DOF::U1>(const shared_ptr<Node>& t_node) { return t_node->get_current_displacement().head(1); }

template<> inline vec get_current_velocity<DOF::U1>(const shared_ptr<Node>& t_node) { return t_node->get_current_velocity().head(1); }

template<> inline vec get_current_acceleration<DOF::U1>(const shared_ptr<Node>& t_node) { return t_node->get_current_acceleration().head(1); }

template<> inline vec get_trial_displacement<DOF::U1>(const shared_ptr<Node>& t_node) { return t_node->get_trial_displacement().head(1); }

template<> inline vec get_trial_velocity<DOF::U1>(const shared_ptr<Node>& t_node) { return t_node->get_trial_velocity().head(1); }

template<> inline vec get_trial_acceleration<DOF::U1>(const shared_ptr<Node>& t_node) { return t_node->get_trial_acceleration().head(1); }

template<> inline vec get_incre_displacement<DOF::U1>(const shared_ptr<Node>& t_node) { return t_node->get_incre_displacement().head(1); }

template<> inline vec get_incre_velocity<DOF::U1>(const shared_ptr<Node>& t_node) { return t_node->get_incre_velocity().head(1); }

template<> inline vec get_incre_acceleration<DOF::U1>(const shared_ptr<Node>& t_node) { return t_node->get_incre_acceleration().head(1); }

template<> inline vec get_current_position<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) {
    const auto& t_coor = t_node->get_coordinate();
    vec t_vec = t_node->get_current_displacement().head(2);
    for(auto I = 0llu; I < std::min(2llu, t_coor.n_elem); ++I) t_vec(I) += t_coor(I);
    return t_vec;
}

template<> inline vec get_trial_position<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) {
    const auto& t_coor = t_node->get_coordinate();
    vec t_vec = t_node->get_trial_displacement().head(2);
    for(auto I = 0llu; I < std::min(2llu, t_coor.n_elem); ++I) t_vec(I) += t_coor(I);
    return t_vec;
}

template<> inline vec get_current_displacement<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) { return t_node->get_current_displacement().head(2); }

template<> inline vec get_current_velocity<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) { return t_node->get_current_velocity().head(2); }

template<> inline vec get_current_acceleration<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) { return t_node->get_current_acceleration().head(2); }

template<> inline vec get_trial_displacement<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) { return t_node->get_trial_displacement().head(2); }

template<> inline vec get_trial_velocity<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) { return t_node->get_trial_velocity().head(2); }

template<> inline vec get_trial_acceleration<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) { return t_node->get_trial_acceleration().head(2); }

template<> inline vec get_incre_displacement<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) { return t_node->get_incre_displacement().head(2); }

template<> inline vec get_incre_velocity<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) { return t_node->get_incre_velocity().head(2); }

template<> inline vec get_incre_acceleration<DOF::U1, DOF::U2>(const shared_ptr<Node>& t_node) { return t_node->get_incre_acceleration().head(2); }

template<> inline vec get_current_position<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) {
    const auto& t_coor = t_node->get_coordinate();
    vec t_vec = t_node->get_current_displacement().head(3);
    for(auto I = 0llu; I < std::min(3llu, t_coor.n_elem); ++I) t_vec(I) += t_coor(I);
    return t_vec;
}

template<> inline vec get_trial_position<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) {
    const auto& t_coor = t_node->get_coordinate();
    vec t_vec = t_node->get_trial_displacement().head(3);
    for(auto I = 0llu; I < std::min(3llu, t_coor.n_elem); ++I) t_vec(I) += t_coor(I);
    return t_vec;
}

template<> inline vec get_current_displacement<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) { return t_node->get_current_displacement().head(3); }

template<> inline vec get_current_velocity<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) { return t_node->get_current_velocity().head(3); }

template<> inline vec get_current_acceleration<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) { return t_node->get_current_acceleration().head(3); }

template<> inline vec get_trial_displacement<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) { return t_node->get_trial_displacement().head(3); }

template<> inline vec get_trial_velocity<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) { return t_node->get_trial_velocity().head(3); }

template<> inline vec get_trial_acceleration<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) { return t_node->get_trial_acceleration().head(3); }

template<> inline vec get_incre_displacement<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) { return t_node->get_incre_displacement().head(3); }

template<> inline vec get_incre_velocity<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) { return t_node->get_incre_velocity().head(3); }

template<> inline vec get_incre_acceleration<DOF::U1, DOF::U2, DOF::U3>(const shared_ptr<Node>& t_node) { return t_node->get_incre_acceleration().head(3); }

#endif
