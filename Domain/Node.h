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
 * @class Node
 * @brief The Node class holds the number of DoFs, coordinate, displacement,
 * velocity and acceleration.
 *
 * The current/committed, incremental and trial status of displacement, velocity
 * and acceleration are stored. These variables will be the communication
 * bridge(s) between Domain, Workshop and Element objects. That is, Element
 * objects do not directly get information from the Workshop. Instead, Workshop
 * passes information to Node objects through the Domain, Element objects
 * acquire new status from associated Node objects only. In this manner, the
 * relationship between those modules remains simple.
 *
 * @author tlc
 * @date 22/07/2017
 * @version 0.2.0
 * @file Node.h
 * @addtogroup Domain
 * @{
 */

#ifndef NODE_H
#define NODE_H

#include <Domain/Tag.h>

class DomainBase;
enum class OutputType;

struct NodeData {
    unsigned num_dof{0u}; // number of DoFs

    vec coordinate; // coordinates of the node

    uvec original_dof;  // original indices
    uvec reordered_dof; // renumbered indices

    vec current_resistance;       // current resistance
    vec current_damping_force;    // current damping force
    vec current_nonviscous_force; // current damping force
    vec current_inertial_force;   // current inertial force
    vec current_displacement;     // current displacement
    vec current_velocity;         // current velocity
    vec current_acceleration;     // current acceleration

    vec incre_resistance;       // incremental resistance
    vec incre_damping_force;    // incremental damping force
    vec incre_nonviscous_force; // incremental damping force
    vec incre_inertial_force;   // incremental inertial force
    vec incre_displacement;     // incremental displacement
    vec incre_velocity;         // incremental velocity
    vec incre_acceleration;     // incremental acceleration

    vec trial_resistance;       // trial resistance
    vec trial_damping_force;    // trial damping force
    vec trial_nonviscous_force; // trial damping force
    vec trial_inertial_force;   // trial inertial force
    vec trial_displacement;     // trial displacement
    vec trial_velocity;         // trial velocity
    vec trial_acceleration;     // trial acceleration
};

class Node final : protected NodeData, public UniqueTag {
public:
    enum class DOF : std::uint8_t {
        NONE,
        U1,          // displacement in x direction
        U2,          // displacement in y direction
        U3,          // displacement in z direction
        UR1,         // rotation in x direction
        UR2,         // rotation in y direction
        UR3,         // rotation in z direction
        FU1,         // fluid displacement
        FU2,         // fluid displacement
        FU3,         // fluid displacement
        FUR1,        // fluid rotation
        FUR2,        // fluid rotation
        FUR3,        // fluid rotation
        RADIAL,      // radius direction in axis-symmetric problem
        AXIAL,       // axial direction in axis-symmetric problem or single section elements
        RS,          // strong axis rotation
        RW,          // weak axis rotation
        DAMAGE,      // damage
        PRESSURE,    // pressure
        TEMPERATURE, // temperature
        WARP         // warping
    };

private:
    std::mutex node_mutex;

    std::vector<DOF> dof_identifier;

public:
    Node(unsigned, vec&&);

    void initialize(const shared_ptr<DomainBase>&);

    void deinitialize();

    void ensure_dof(unsigned, const std::vector<DOF>&);
    [[nodiscard]] bool validate_dof(const std::vector<DOF>&) const;
    [[nodiscard]] std::vector<uword> get_dof(const std::vector<DOF>&) const;

    void set_original_dof(unsigned&);
    [[nodiscard]] const uvec& get_original_dof() const;

    void set_reordered_dof(const uvec&);
    [[nodiscard]] const uvec& get_reordered_dof() const;

    [[nodiscard]] const vec& get_coordinate() const;

    [[nodiscard]] const vec& get_current_resistance() const;
    [[nodiscard]] const vec& get_current_damping_force() const;
    [[nodiscard]] const vec& get_current_nonviscous_force() const;
    [[nodiscard]] const vec& get_current_inertial_force() const;
    [[nodiscard]] const vec& get_current_displacement() const;
    [[nodiscard]] const vec& get_current_velocity() const;
    [[nodiscard]] const vec& get_current_acceleration() const;

    [[nodiscard]] const vec& get_incre_resistance() const;
    [[nodiscard]] const vec& get_incre_damping_force() const;
    [[nodiscard]] const vec& get_incre_nonviscous_force() const;
    [[nodiscard]] const vec& get_incre_inertial_force() const;
    [[nodiscard]] const vec& get_incre_displacement() const;
    [[nodiscard]] const vec& get_incre_velocity() const;
    [[nodiscard]] const vec& get_incre_acceleration() const;

    [[nodiscard]] const vec& get_trial_resistance() const;
    [[nodiscard]] const vec& get_trial_damping_force() const;
    [[nodiscard]] const vec& get_trial_nonviscous_force() const;
    [[nodiscard]] const vec& get_trial_inertial_force() const;
    [[nodiscard]] const vec& get_trial_displacement() const;
    [[nodiscard]] const vec& get_trial_velocity() const;
    [[nodiscard]] const vec& get_trial_acceleration() const;

    [[nodiscard]] vec get_current_resistance(unsigned) const;
    [[nodiscard]] vec get_current_damping_force(unsigned) const;
    [[nodiscard]] vec get_current_nonviscous_force(unsigned) const;
    [[nodiscard]] vec get_current_inertial_force(unsigned) const;
    [[nodiscard]] vec get_current_displacement(unsigned) const;
    [[nodiscard]] vec get_current_velocity(unsigned) const;
    [[nodiscard]] vec get_current_acceleration(unsigned) const;

    [[nodiscard]] vec get_incre_resistance(unsigned) const;
    [[nodiscard]] vec get_incre_damping_force(unsigned) const;
    [[nodiscard]] vec get_incre_nonviscous_force(unsigned) const;
    [[nodiscard]] vec get_incre_inertial_force(unsigned) const;
    [[nodiscard]] vec get_incre_displacement(unsigned) const;
    [[nodiscard]] vec get_incre_velocity(unsigned) const;
    [[nodiscard]] vec get_incre_acceleration(unsigned) const;

    [[nodiscard]] vec get_trial_resistance(unsigned) const;
    [[nodiscard]] vec get_trial_damping_force(unsigned) const;
    [[nodiscard]] vec get_trial_nonviscous_force(unsigned) const;
    [[nodiscard]] vec get_trial_inertial_force(unsigned) const;
    [[nodiscard]] vec get_trial_displacement(unsigned) const;
    [[nodiscard]] vec get_trial_velocity(unsigned) const;
    [[nodiscard]] vec get_trial_acceleration(unsigned) const;

    [[nodiscard]] vec initial_position(unsigned) const;
    [[nodiscard]] vec current_position(unsigned) const;
    [[nodiscard]] vec trial_position(unsigned) const;

    const vec& update_current_resistance(vec&&);
    const vec& update_current_damping_force(vec&&);
    const vec& update_current_nonviscous_force(vec&&);
    const vec& update_current_inertial_force(vec&&);
    const vec& update_current_displacement(vec&&);
    const vec& update_current_velocity(vec&&);
    const vec& update_current_acceleration(vec&&);

    const vec& update_incre_resistance(vec&&);
    const vec& update_incre_damping_force(vec&&);
    const vec& update_incre_nonviscous_force(vec&&);
    const vec& update_incre_inertial_force(vec&&);
    const vec& update_incre_displacement(vec&&);
    const vec& update_incre_velocity(vec&&);
    const vec& update_incre_acceleration(vec&&);

    const vec& update_trial_resistance(vec&&);
    const vec& update_trial_damping_force(vec&&);
    const vec& update_trial_nonviscous_force(vec&&);
    const vec& update_trial_inertial_force(vec&&);
    const vec& update_trial_displacement(vec&&);
    const vec& update_trial_velocity(vec&&);
    const vec& update_trial_acceleration(vec&&);

    void update_current_status(const vec&);
    void update_current_status(const vec&, const vec&);
    void update_current_status(const vec&, const vec&, const vec&);

    void update_incre_status(const vec&);
    void update_incre_status(const vec&, const vec&);
    void update_incre_status(const vec&, const vec&, const vec&);

    void update_trial_status(const vec&);
    void update_trial_status(const vec&, const vec&);
    void update_trial_status(const vec&, const vec&, const vec&);

    void commit_status();
    void reset_status();
    void clear_status();

    [[nodiscard]] std::vector<vec> record(OutputType) const;

    void print() override;
};

namespace suanpan {
#if defined(__GNUC__) && (__GNUC__ < 12)
    inline
#else
    constexpr inline
#endif
        std::vector<Node::DOF>
        translational(const unsigned dimension) {
        return 2u == dimension ? std::vector{Node::DOF::U1, Node::DOF::U2} : std::vector{Node::DOF::U1, Node::DOF::U2, Node::DOF::U3};
    }

#if defined(__GNUC__) && (__GNUC__ < 12)
    inline
#else
    constexpr inline
#endif
        std::vector<Node::DOF>
        mechanical(const unsigned dimension) {
        return 2u == dimension ? std::vector{Node::DOF::U1, Node::DOF::U2, Node::DOF::UR3} : std::vector{Node::DOF::U1, Node::DOF::U2, Node::DOF::U3, Node::DOF::UR1, Node::DOF::UR2, Node::DOF::UR3};
    }
} // namespace suanpan

#endif

//! @}
