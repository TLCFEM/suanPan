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
 * @class Element
 * @brief A Element class.
 * @author tlc
 * @date 12/01/2020
 * @version 0.3.0
 * @file Element.h
 * @addtogroup Element
 * @{
 */

#ifndef ELEMENT_H
#define ELEMENT_H

#include <Element/Distributed.h>
#include <Element/ElementBase.h>

enum class MaterialType : unsigned;
enum class SectionType : unsigned;

//struct StateElement {
//    mat mass{};
//    mat damping{};
//    mat stiffness{};
//    mat geometry{};
//
//    vec resistance{};
//    vec damping_force{};
//    vec inertial_force{};
//
//    vec body_force{};
//    vec traction{};
//};

struct DataElement {
    const uvec node_encoding; // node encoding
    const uvec material_tag;  // material tags
    const uvec section_tag;   // section tags

    const bool nlgeom = false; // nonlinear geometry switch

    bool update_mass = true;       // flag to indicate if update matrix
    bool update_viscous = true;    // flag to indicate if update matrix
    bool update_nonviscous = true; // flag to indicate if update matrix
    bool update_stiffness = true;  // flag to indicate if update matrix
    bool update_geometry = true;   // flag to indicate if update matrix

    bool modify_mass = true;       // flag to indicate if modify matrix
    bool modify_viscous = true;    // flag to indicate if modify matrix
    bool modify_nonviscous = true; // flag to indicate if modify matrix

    uvec dof_encoding{}; // DoF encoding vector

    mat initial_mass{};       // mass matrix
    mat initial_viscous{};    // viscous damping matrix
    mat initial_nonviscous{}; // nonviscous damping matrix
    mat initial_stiffness{};  // stiffness matrix
    mat initial_geometry{};   // geometry matrix

    mat trial_mass{};       // mass matrix
    mat trial_viscous{};    // viscous damping matrix
    mat trial_nonviscous{}; // nonviscous damping matrix
    mat trial_stiffness{};  // stiffness matrix
    mat trial_geometry{};   // geometry matrix

    mat current_mass{};       // mass matrix
    mat current_viscous{};    // viscous damping matrix
    mat current_nonviscous{}; // nonviscous damping matrix
    mat current_stiffness{};  // stiffness matrix
    mat current_geometry{};   // geometry matrix

    vec trial_resistance{};            // resistance vector
    vec current_resistance{};          // resistance vector
    vec trial_viscous_force{};         // viscous damping force
    vec current_viscous_force{};       // viscous damping force
    cx_mat trial_nonviscous_force{};   // nonviscous damping force
    cx_mat current_nonviscous_force{}; // nonviscous damping force
    vec trial_inertial_force{};        // inertial force
    vec current_inertial_force{};      // inertial force

    vec trial_body_force{};
    vec current_body_force{};
    vec trial_traction{};
    vec current_traction{};

    mat body_force{};
    mat traction{};

    mat stiffness_container{}; // universal container to accommodate matrix based on modified stiffness matrix
    mat mass_container{};      // universal container to accommodate matrix based on modified mass matrix

    double strain_energy = 0.;
    double kinetic_energy = 0.;
    double viscous_energy = 0.;
    double nonviscous_energy = 0.;
    double complementary_energy = 0.;
    vec momentum{};

    const double characteristic_length = 1.;
};

class Element : protected DataElement, public ElementBase, public Distributed {
    const unsigned num_node;                      // number of nodes
    const unsigned num_dof;                       // number of DoFs
    const unsigned num_size = num_dof * num_node; // number of size

    const bool initialized = false;
    const bool symmetric = false;
    const bool use_group = false;
    const unsigned use_other = 0;

    const MaterialType material_type;
    const SectionType section_type;

    const std::vector<DOF> dof_identifier;

    std::vector<MappingDOF> dof_mapping;

    friend void ConstantMass(DataElement*);
    friend void ConstantDamping(DataElement*);
    friend void ConstantStiffness(DataElement*);
    friend void ConstantGeometry(DataElement*);

    void update_strain_energy() override;
    void update_kinetic_energy() override;
    void update_viscous_energy() override;
    void update_nonviscous_energy() override;
    void update_complementary_energy() override;
    void update_momentum() override;

protected:
    std::vector<weak_ptr<Node>> node_ptr; // node pointers

    [[nodiscard]] mat get_coordinate(unsigned) const override;

    [[nodiscard]] vec get_node_incre_resistance() const override;
    [[nodiscard]] vec get_node_trial_resistance() const override;
    [[nodiscard]] vec get_node_current_resistance() const override;

    [[nodiscard]] std::vector<shared_ptr<Material>> get_material(const shared_ptr<DomainBase>&) const override;
    [[nodiscard]] std::vector<shared_ptr<Section>> get_section(const shared_ptr<DomainBase>&) const override;

public:
    Element(
        unsigned,          // tag
        unsigned,          // number of nodes
        unsigned,          // number of dofs
        uvec&&,            // node encoding
        std::vector<DOF>&& // dof identifier
    );
    Element(
        unsigned,          // tag
        unsigned,          // number of nodes
        unsigned,          // number of dofs
        uvec&&,            // node encoding
        uvec&&,            // material tags
        bool,              // nonlinear geometry switch
        MaterialType,      // material type for internal check
        std::vector<DOF>&& // dof identifier
    );
    Element(
        unsigned,          // tag
        unsigned,          // number of nodes
        unsigned,          // number of dofs
        uvec&&,            // node encoding
        uvec&&,            // section tags
        bool,              // nonlinear geometry switch
        SectionType,       // section type for internal check
        std::vector<DOF>&& // dof identifier
    );
    Element(
        unsigned, // tag
        unsigned, // number of dofs
        uvec&&    // group encoding
    );
    Element(
        unsigned, // tag
        unsigned, // number of dofs
        unsigned, // other element tag
        unsigned  // node tag
    );
    Element(const Element&) = delete;            // copy forbidden
    Element(Element&&) = delete;                 // move forbidden
    Element& operator=(const Element&) = delete; // assign forbidden
    Element& operator=(Element&&) = delete;      // assign forbidden

    ~Element() override = default;

    int initialize_base(const shared_ptr<DomainBase>&) final;

    void set_initialized(bool) const override;
    void set_symmetric(bool) const override;
    [[nodiscard]] bool is_initialized() const override;
    [[nodiscard]] bool is_symmetric() const override;
    [[nodiscard]] bool is_nlgeom() const override;

    void update_dof_encoding() override;

    [[nodiscard]] bool if_update_mass() const override;
    [[nodiscard]] bool if_update_viscous() const override;
    [[nodiscard]] bool if_update_nonviscous() const override;
    [[nodiscard]] bool if_update_stiffness() const override;
    [[nodiscard]] bool if_update_geometry() const override;

    [[nodiscard]] bool allow_modify_mass() const override;
    [[nodiscard]] bool allow_modify_viscous() const override;
    [[nodiscard]] bool allow_modify_nonviscous() const override;

    [[nodiscard]] const uvec& get_dof_encoding() const override;
    [[nodiscard]] const uvec& get_node_encoding() const override;

    [[nodiscard]] const std::vector<MappingDOF>& get_dof_mapping() const override;

    [[nodiscard]] const uvec& get_material_tag() const override;
    [[nodiscard]] const uvec& get_section_tag() const override;

    [[nodiscard]] unsigned get_dof_number() const override;
    [[nodiscard]] unsigned get_node_number() const override;
    [[nodiscard]] unsigned get_total_number() const override;

    void clear_node_ptr() override;
    [[nodiscard]] const std::vector<weak_ptr<Node>>& get_node_ptr() const override;

    [[nodiscard]] vec get_incre_displacement() const override;
    [[nodiscard]] vec get_incre_velocity() const override;
    [[nodiscard]] vec get_incre_acceleration() const override;
    [[nodiscard]] vec get_trial_displacement() const override;
    [[nodiscard]] vec get_trial_velocity() const override;
    [[nodiscard]] vec get_trial_acceleration() const override;
    [[nodiscard]] vec get_current_displacement() const override;
    [[nodiscard]] vec get_current_velocity() const override;
    [[nodiscard]] vec get_current_acceleration() const override;

    [[nodiscard]] const vec& get_trial_resistance() const override;
    [[nodiscard]] const vec& get_current_resistance() const override;
    [[nodiscard]] const vec& get_trial_damping_force() const override;
    [[nodiscard]] const vec& get_current_damping_force() const override;
    [[nodiscard]] const cx_mat& get_trial_nonviscous_force() const override;
    [[nodiscard]] const cx_mat& get_current_nonviscous_force() const override;
    [[nodiscard]] const vec& get_trial_inertial_force() override;
    [[nodiscard]] const vec& get_current_inertial_force() override;

    [[nodiscard]] const vec& get_trial_body_force() const override;
    [[nodiscard]] const vec& get_current_body_force() const override;
    [[nodiscard]] const vec& get_trial_traction() const override;
    [[nodiscard]] const vec& get_current_traction() const override;

    [[nodiscard]] const mat& get_trial_mass() const override;
    [[nodiscard]] const mat& get_trial_viscous() const override;
    [[nodiscard]] const mat& get_trial_nonviscous() const override;
    [[nodiscard]] const mat& get_trial_stiffness() const override;
    [[nodiscard]] const mat& get_trial_geometry() const override;
    [[nodiscard]] const mat& get_trial_secant() const override;

    [[nodiscard]] const mat& get_current_mass() const override;
    [[nodiscard]] const mat& get_current_viscous() const override;
    [[nodiscard]] const mat& get_current_nonviscous() const override;
    [[nodiscard]] const mat& get_current_stiffness() const override;
    [[nodiscard]] const mat& get_current_geometry() const override;
    [[nodiscard]] const mat& get_current_secant() const override;

    [[nodiscard]] const mat& get_initial_mass() const override;
    [[nodiscard]] const mat& get_initial_viscous() const override;
    [[nodiscard]] const mat& get_initial_nonviscous() const override;
    [[nodiscard]] const mat& get_initial_stiffness() const override;
    [[nodiscard]] const mat& get_initial_geometry() const override;
    [[nodiscard]] const mat& get_initial_secant() const override;

    [[nodiscard]] const mat& get_mass_container() const override;
    [[nodiscard]] const mat& get_stiffness_container() const override;

    int clear_status() override = 0;
    int commit_status() override = 0;
    int reset_status() override = 0;

    const vec& update_body_force(const vec&) override;
    const vec& update_traction(const vec&) override;

    std::vector<vec> record(OutputType) override;

    [[nodiscard]] double get_strain_energy() const override;
    [[nodiscard]] double get_complementary_energy() const override;
    [[nodiscard]] double get_kinetic_energy() const override;
    [[nodiscard]] double get_viscous_energy() const override;
    [[nodiscard]] double get_nonviscous_energy() const override;
    [[nodiscard]] const vec& get_momentum() const override;
    [[nodiscard]] double get_momentum_component(DOF) const override;

    [[nodiscard]] double get_characteristic_length() const override;

    [[nodiscard]] mat compute_shape_function(const mat&, unsigned) const override;
};

std::vector<vec>& append_to(std::vector<vec>&, std::vector<vec>&&);

#endif

//! @}
