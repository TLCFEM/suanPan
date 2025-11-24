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

#include "Node.h"

#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>
#include <Toolbox/utility.h>

Node::Node(const unsigned T, vec&& C)
    : UniqueTag(T) { coordinate = std::move(C); }

/**
 * \brief This method should be called after elements are set.
 * Elements will set the minimum number of DoFs for all related nodes.
 */
void Node::initialize(const shared_ptr<DomainBase>& D) {
    if(!is_active()) return;

    if(0u == num_dof) {
        suanpan_debug("Node {} disabled as it is not used.\n", get_tag());
        D->disable_node(get_tag());
        return;
    }

    original_dof.reset();
    reordered_dof.reset();
    dof_identifier.clear();

    current_displacement.resize(num_dof);
    current_velocity.resize(num_dof);
    current_acceleration.resize(num_dof);

    incre_displacement.resize(num_dof);
    incre_velocity.resize(num_dof);
    incre_acceleration.resize(num_dof);

    trial_displacement.resize(num_dof);
    trial_velocity.resize(num_dof);
    trial_acceleration.resize(num_dof);
}

void Node::deinitialize() { num_dof = 0u; }

void Node::ensure_dof_number(const unsigned D) {
    std::scoped_lock node_lock{node_mutex};

    if(num_dof >= D) return;

    num_dof = D;
}

void Node::set_dof_identifier(const std::vector<DOF>& D) {
    std::scoped_lock node_lock{node_mutex};

    if(dof_identifier.empty()) dof_identifier.resize(num_dof, DOF::NONE);

    for(size_t I = 0; I < D.size(); ++I) {
        if(DOF::NONE == D[I]) continue;
        if(DOF::NONE != dof_identifier[I] && D[I] != dof_identifier[I])
            suanpan_warning("Inconsistent DoF assignment for node {} detected.\n", get_tag());
        dof_identifier[I] = D[I];
    }
}

const std::vector<Node::DOF>& Node::get_dof_identifier() const { return dof_identifier; }

void Node::set_original_dof(unsigned& F) {
    original_dof.set_size(num_dof);
    for(auto& I : original_dof) I = F++;
}

const uvec& Node::get_original_dof() const { return original_dof; }

void Node::set_reordered_dof(const uvec& R) { reordered_dof = R(original_dof); }

const uvec& Node::get_reordered_dof() const { return reordered_dof.is_empty() ? original_dof : reordered_dof; }

const vec& Node::get_coordinate() const { return coordinate; }

const vec& Node::get_current_resistance() const { return current_resistance; }

const vec& Node::get_current_damping_force() const { return current_damping_force; }

const vec& Node::get_current_nonviscous_force() const { return current_nonviscous_force; }

const vec& Node::get_current_inertial_force() const { return current_inertial_force; }

const vec& Node::get_current_displacement() const { return current_displacement; }

const vec& Node::get_current_velocity() const { return current_velocity; }

const vec& Node::get_current_acceleration() const { return current_acceleration; }

const vec& Node::get_incre_resistance() const { return incre_resistance; }

const vec& Node::get_incre_damping_force() const { return incre_damping_force; }

const vec& Node::get_incre_nonviscous_force() const { return incre_nonviscous_force; }

const vec& Node::get_incre_inertial_force() const { return incre_inertial_force; }

const vec& Node::get_incre_displacement() const { return incre_displacement; }

const vec& Node::get_incre_velocity() const { return incre_velocity; }

const vec& Node::get_incre_acceleration() const { return incre_acceleration; }

const vec& Node::get_trial_resistance() const { return trial_resistance; }

const vec& Node::get_trial_damping_force() const { return trial_damping_force; }

const vec& Node::get_trial_nonviscous_force() const { return trial_nonviscous_force; }

const vec& Node::get_trial_inertial_force() const { return trial_inertial_force; }

const vec& Node::get_trial_displacement() const { return trial_displacement; }

const vec& Node::get_trial_velocity() const { return trial_velocity; }

const vec& Node::get_trial_acceleration() const { return trial_acceleration; }

const vec& Node::update_current_resistance(vec&& in) {
    trial_resistance = current_resistance = std::move(in);
    incre_resistance.zeros(current_resistance.size());
    return current_resistance;
}

const vec& Node::update_current_damping_force(vec&& in) {
    trial_damping_force = current_damping_force = std::move(in);
    incre_damping_force.zeros(current_damping_force.size());
    return current_damping_force;
}

const vec& Node::update_current_nonviscous_force(vec&& in) {
    trial_nonviscous_force = current_nonviscous_force = std::move(in);
    incre_nonviscous_force.zeros(current_nonviscous_force.size());
    return current_nonviscous_force;
}

const vec& Node::update_current_inertial_force(vec&& in) {
    trial_inertial_force = current_inertial_force = std::move(in);
    incre_inertial_force.zeros(current_inertial_force.size());
    return current_inertial_force;
}

const vec& Node::update_current_displacement(vec&& in) {
    trial_displacement = current_displacement = std::move(in);
    incre_displacement.zeros(current_displacement.size());
    return current_displacement;
}

const vec& Node::update_current_velocity(vec&& in) {
    trial_velocity = current_velocity = std::move(in);
    incre_velocity.zeros(current_velocity.size());
    return current_velocity;
}

const vec& Node::update_current_acceleration(vec&& in) {
    trial_acceleration = current_acceleration = std::move(in);
    incre_acceleration.zeros(current_acceleration.size());
    return current_acceleration;
}

const vec& Node::update_incre_resistance(vec&& in) {
    incre_resistance = std::move(in);
    trial_resistance = current_resistance.resize(incre_resistance.size()) + incre_resistance;
    return incre_resistance;
}

const vec& Node::update_incre_damping_force(vec&& in) {
    incre_damping_force = std::move(in);
    trial_damping_force = current_damping_force.resize(incre_damping_force.size()) + incre_damping_force;
    return incre_damping_force;
}

const vec& Node::update_incre_nonviscous_force(vec&& in) {
    incre_nonviscous_force = std::move(in);
    trial_nonviscous_force = current_nonviscous_force.resize(incre_nonviscous_force.size()) + incre_nonviscous_force;
    return incre_nonviscous_force;
}

const vec& Node::update_incre_inertial_force(vec&& in) {
    incre_inertial_force = std::move(in);
    trial_inertial_force = current_inertial_force.resize(incre_inertial_force.size()) + incre_inertial_force;
    return incre_inertial_force;
}

const vec& Node::update_incre_displacement(vec&& in) {
    incre_displacement = std::move(in);
    trial_displacement = current_displacement.resize(incre_displacement.size()) + incre_displacement;
    return incre_displacement;
}

const vec& Node::update_incre_velocity(vec&& in) {
    incre_velocity = std::move(in);
    trial_velocity = current_velocity.resize(incre_velocity.size()) + incre_velocity;
    return incre_velocity;
}

const vec& Node::update_incre_acceleration(vec&& in) {
    incre_acceleration = std::move(in);
    trial_acceleration = current_acceleration.resize(incre_acceleration.size()) + incre_acceleration;
    return incre_acceleration;
}

const vec& Node::update_trial_resistance(vec&& in) {
    trial_resistance = std::move(in);
    incre_resistance = trial_resistance - current_resistance.resize(trial_resistance.size());
    return trial_resistance;
}

const vec& Node::update_trial_damping_force(vec&& in) {
    trial_damping_force = std::move(in);
    incre_damping_force = trial_damping_force - current_damping_force.resize(trial_damping_force.size());
    return trial_damping_force;
}

const vec& Node::update_trial_nonviscous_force(vec&& in) {
    trial_nonviscous_force = std::move(in);
    incre_nonviscous_force = trial_nonviscous_force - current_nonviscous_force.resize(trial_nonviscous_force.size());
    return trial_nonviscous_force;
}

const vec& Node::update_trial_inertial_force(vec&& in) {
    trial_inertial_force = std::move(in);
    incre_inertial_force = trial_inertial_force - current_inertial_force.resize(trial_inertial_force.size());
    return trial_inertial_force;
}

const vec& Node::update_trial_displacement(vec&& in) {
    trial_displacement = std::move(in);
    incre_displacement = trial_displacement - current_displacement.resize(trial_displacement.size());
    return trial_displacement;
}

const vec& Node::update_trial_velocity(vec&& in) {
    trial_velocity = std::move(in);
    incre_velocity = trial_velocity - current_velocity.resize(trial_velocity.size());
    return trial_velocity;
}

const vec& Node::update_trial_acceleration(vec&& in) {
    trial_acceleration = std::move(in);
    incre_acceleration = trial_acceleration - current_acceleration.resize(trial_acceleration.size());
    return trial_acceleration;
}

void Node::update_current_status(const vec& D) { update_current_displacement(D(reordered_dof)); }

void Node::update_current_status(const vec& D, const vec& V) {
    update_current_velocity(V(reordered_dof));
    update_current_status(D);
}

void Node::update_current_status(const vec& D, const vec& V, const vec& A) {
    update_current_acceleration(A(reordered_dof));
    update_current_status(D, V);
}

void Node::update_incre_status(const vec& D) { update_incre_displacement(D(reordered_dof)); }

void Node::update_incre_status(const vec& D, const vec& V) {
    update_incre_velocity(V(reordered_dof));
    update_incre_status(D);
}

void Node::update_incre_status(const vec& D, const vec& V, const vec& A) {
    update_incre_acceleration(A(reordered_dof));
    update_incre_status(D, V);
}

void Node::update_trial_status(const vec& D) { update_trial_displacement(D(reordered_dof)); }

void Node::update_trial_status(const vec& D, const vec& V) {
    update_trial_velocity(V(reordered_dof));
    update_trial_status(D);
}

void Node::update_trial_status(const vec& D, const vec& V, const vec& A) {
    update_trial_acceleration(A(reordered_dof));
    update_trial_status(D, V);
}

void Node::commit_status() {
    if(!trial_resistance.is_empty()) {
        current_resistance = trial_resistance;
        incre_resistance.zeros();
    }
    if(!trial_damping_force.is_empty()) {
        current_damping_force = trial_damping_force;
        incre_damping_force.zeros();
    }
    if(!trial_nonviscous_force.is_empty()) {
        current_nonviscous_force = trial_nonviscous_force;
        incre_nonviscous_force.zeros();
    }
    if(!trial_inertial_force.is_empty()) {
        current_inertial_force = trial_inertial_force;
        incre_inertial_force.zeros();
    }
    if(!trial_displacement.is_empty()) {
        current_displacement = trial_displacement;
        incre_displacement.zeros();
    }
    if(!trial_velocity.is_empty()) {
        current_velocity = trial_velocity;
        incre_displacement.zeros();
    }
    if(!trial_acceleration.is_empty()) {
        current_acceleration = trial_acceleration;
        incre_acceleration.zeros();
    }
}

void Node::reset_status() {
    if(!current_resistance.is_empty()) {
        trial_resistance = current_resistance;
        incre_resistance.zeros();
    }
    if(!current_damping_force.is_empty()) {
        trial_damping_force = current_damping_force;
        incre_damping_force.zeros();
    }
    if(!current_nonviscous_force.is_empty()) {
        trial_nonviscous_force = current_nonviscous_force;
        incre_nonviscous_force.zeros();
    }
    if(!current_inertial_force.is_empty()) {
        trial_inertial_force = current_inertial_force;
        incre_inertial_force.zeros();
    }
    if(!current_displacement.is_empty()) {
        trial_displacement = current_displacement;
        incre_displacement.zeros();
    }
    if(!current_velocity.is_empty()) {
        trial_velocity = current_velocity;
        incre_velocity.zeros();
    }
    if(!current_acceleration.is_empty()) {
        trial_acceleration = current_acceleration;
        incre_acceleration.zeros();
    }
}

void Node::clear_status() {
    if(!current_resistance.is_empty()) {
        current_resistance.zeros();
        incre_resistance.zeros();
        trial_resistance.zeros();
    }
    if(!current_damping_force.is_empty()) {
        current_damping_force.zeros();
        incre_damping_force.zeros();
        trial_damping_force.zeros();
    }
    if(!current_nonviscous_force.is_empty()) {
        current_nonviscous_force.zeros();
        incre_nonviscous_force.zeros();
        trial_nonviscous_force.zeros();
    }
    if(!current_inertial_force.is_empty()) {
        current_inertial_force.zeros();
        incre_inertial_force.zeros();
        trial_inertial_force.zeros();
    }
    if(!current_displacement.is_empty()) {
        current_displacement.zeros();
        incre_displacement.zeros();
        trial_displacement.zeros();
    }
    if(!current_velocity.is_empty()) {
        current_velocity.zeros();
        incre_velocity.zeros();
        trial_velocity.zeros();
    }
    if(!current_acceleration.is_empty()) {
        current_acceleration.zeros();
        incre_acceleration.zeros();
        trial_acceleration.zeros();
    }
}

std::vector<vec> Node::record(const OutputType L) const {
    std::vector<vec> data;

    const auto ensure = [&](const vec& in) { return in.empty() ? zeros(num_dof) : in; };

    if(L == OutputType::RF) data.emplace_back(ensure(current_resistance));
    else if(L == OutputType::DF) data.emplace_back(ensure(current_damping_force));
    else if(L == OutputType::IF) data.emplace_back(ensure(current_inertial_force));
    else if(L == OutputType::U) data.emplace_back(ensure(current_displacement));
    else if(L == OutputType::V) data.emplace_back(ensure(current_velocity));
    else if(L == OutputType::A) data.emplace_back(ensure(current_acceleration));

    return data;
}

void Node::print() {
    suanpan_info("Node {}{}\n", get_tag(), is_active() ? ":" : " is currently inactive.");
    suanpan_info("Coordinate:", coordinate);
    suanpan_info("Displacement:", current_displacement);
    suanpan_info("Resistance:", current_resistance);
    if(!suanpan::approx_equal(accu(current_velocity), 0.))
        suanpan_info("Velocity:", current_velocity);
    if(!suanpan::approx_equal(accu(current_acceleration), 0.))
        suanpan_info("Acceleration:", current_acceleration);
}
