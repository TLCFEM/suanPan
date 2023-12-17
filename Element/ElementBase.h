/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class ElementBase
 * @brief A ElementBase class.
 * @author tlc
 * @date 10/01/2020
 * @version 0.1.0
 * @file ElementBase.h
 * @addtogroup Element
 * @{
 */

#ifndef ELEMENTBASE_H
#define ELEMENTBASE_H

#include <Domain/Tag.h>
#include <Element/Visualisation/vtkBase.h>
#include <Element/MappingDOF.h>

class Node;
class DomainBase;
class Material;
class Section;
enum class OutputType;
enum class DOF : unsigned short;

class ElementBase : public Tag, public vtkBase {
    virtual void update_strain_energy() = 0;
    virtual void update_kinetic_energy() = 0;
    virtual void update_viscous_energy() = 0;
    virtual void update_nonviscous_energy() = 0;
    virtual void update_complementary_energy() = 0;
    virtual void update_momentum() = 0;

protected:
    friend mat get_coordinate(const ElementBase*, unsigned);

    friend vec get_incre_displacement(const ElementBase*);
    friend vec get_incre_velocity(const ElementBase*);
    friend vec get_incre_acceleration(const ElementBase*);
    friend vec get_trial_displacement(const ElementBase*);
    friend vec get_trial_velocity(const ElementBase*);
    friend vec get_trial_acceleration(const ElementBase*);
    friend vec get_current_displacement(const ElementBase*);
    friend vec get_current_velocity(const ElementBase*);
    friend vec get_current_acceleration(const ElementBase*);

    [[nodiscard]] virtual mat get_coordinate(unsigned) const = 0;

    [[nodiscard]] virtual vec get_node_incre_resistance() const = 0;
    [[nodiscard]] virtual vec get_node_trial_resistance() const = 0;
    [[nodiscard]] virtual vec get_node_current_resistance() const = 0;

    [[nodiscard]] virtual std::vector<shared_ptr<Material>> get_material(const shared_ptr<DomainBase>&) const = 0;
    [[nodiscard]] virtual std::vector<shared_ptr<Section>> get_section(const shared_ptr<DomainBase>&) const = 0;

public:
    explicit ElementBase(const unsigned T)
        : Tag(T) {}

    ElementBase(const ElementBase&) = delete;            // copy forbidden
    ElementBase(ElementBase&&) = delete;                 // move forbidden
    ElementBase& operator=(const ElementBase&) = delete; // assign forbidden
    ElementBase& operator=(ElementBase&&) = delete;      // assign forbidden

    ~ElementBase() override = default;

    virtual int initialize_base(const shared_ptr<DomainBase>&) = 0;
    virtual int initialize(const shared_ptr<DomainBase>&) = 0;

    virtual void set_initialized(bool) const = 0;
    virtual void set_symmetric(bool) const = 0;
    [[nodiscard]] virtual bool is_initialized() const = 0;
    [[nodiscard]] virtual bool is_symmetric() const = 0;
    [[nodiscard]] virtual bool is_nlgeom() const = 0;

    virtual void update_dof_encoding() = 0;

    [[nodiscard]] virtual bool if_update_mass() const = 0;
    [[nodiscard]] virtual bool if_update_damping() const = 0;
    [[nodiscard]] virtual bool if_update_nonviscous() const = 0;
    [[nodiscard]] virtual bool if_update_stiffness() const = 0;
    [[nodiscard]] virtual bool if_update_geometry() const = 0;

    [[nodiscard]] virtual const uvec& get_dof_encoding() const = 0;
    [[nodiscard]] virtual const uvec& get_node_encoding() const = 0;

    [[nodiscard]] virtual const std::vector<MappingDOF>& get_dof_mapping() const = 0;

    [[nodiscard]] virtual const uvec& get_material_tag() const = 0;
    [[nodiscard]] virtual const uvec& get_section_tag() const = 0;

    [[nodiscard]] virtual unsigned get_dof_number() const = 0;
    [[nodiscard]] virtual unsigned get_node_number() const = 0;
    [[nodiscard]] virtual unsigned get_total_number() const = 0;

    virtual void clear_node_ptr() = 0;
    [[nodiscard]] virtual const std::vector<weak_ptr<Node>>& get_node_ptr() const = 0;

    [[nodiscard]] virtual vec get_incre_displacement() const = 0;
    [[nodiscard]] virtual vec get_incre_velocity() const = 0;
    [[nodiscard]] virtual vec get_incre_acceleration() const = 0;
    [[nodiscard]] virtual vec get_trial_displacement() const = 0;
    [[nodiscard]] virtual vec get_trial_velocity() const = 0;
    [[nodiscard]] virtual vec get_trial_acceleration() const = 0;
    [[nodiscard]] virtual vec get_current_displacement() const = 0;
    [[nodiscard]] virtual vec get_current_velocity() const = 0;
    [[nodiscard]] virtual vec get_current_acceleration() const = 0;

    [[nodiscard]] virtual const vec& get_trial_resistance() const = 0;
    [[nodiscard]] virtual const vec& get_current_resistance() const = 0;
    [[nodiscard]] virtual const vec& get_trial_damping_force() const = 0;
    [[nodiscard]] virtual const vec& get_current_damping_force() const = 0;
    [[nodiscard]] virtual const cx_mat& get_trial_nonviscous_force() const = 0;
    [[nodiscard]] virtual const cx_mat& get_current_nonviscous_force() const = 0;
    [[nodiscard]] virtual const vec& get_trial_inertial_force() = 0;
    [[nodiscard]] virtual const vec& get_current_inertial_force() = 0;

    [[nodiscard]] virtual const vec& get_trial_body_force() const = 0;
    [[nodiscard]] virtual const vec& get_current_body_force() const = 0;
    [[nodiscard]] virtual const vec& get_trial_traction() const = 0;
    [[nodiscard]] virtual const vec& get_current_traction() const = 0;

    [[nodiscard]] virtual const mat& get_trial_mass() const = 0;
    [[nodiscard]] virtual const mat& get_trial_damping() const = 0;
    [[nodiscard]] virtual const mat& get_trial_nonviscous() const = 0;
    [[nodiscard]] virtual const mat& get_trial_stiffness() const = 0;
    [[nodiscard]] virtual const mat& get_trial_geometry() const = 0;
    [[nodiscard]] virtual const mat& get_trial_secant() const = 0;

    [[nodiscard]] virtual const mat& get_current_mass() const = 0;
    [[nodiscard]] virtual const mat& get_current_damping() const = 0;
    [[nodiscard]] virtual const mat& get_current_nonviscous() const = 0;
    [[nodiscard]] virtual const mat& get_current_stiffness() const = 0;
    [[nodiscard]] virtual const mat& get_current_geometry() const = 0;
    [[nodiscard]] virtual const mat& get_current_secant() const = 0;

    [[nodiscard]] virtual const mat& get_initial_mass() const = 0;
    [[nodiscard]] virtual const mat& get_initial_damping() const = 0;
    [[nodiscard]] virtual const mat& get_initial_nonviscous() const = 0;
    [[nodiscard]] virtual const mat& get_initial_stiffness() const = 0;
    [[nodiscard]] virtual const mat& get_initial_geometry() const = 0;
    [[nodiscard]] virtual const mat& get_initial_secant() const = 0;

    [[nodiscard]] virtual const mat& get_mass_container() const = 0;
    [[nodiscard]] virtual const mat& get_stiffness_container() const = 0;

    virtual int update_status() = 0;
    virtual int clear_status() = 0;
    virtual int commit_status() = 0;
    virtual int reset_status() = 0;

    virtual const vec& update_body_force(const vec&) = 0;
    virtual const vec& update_traction(const vec&) = 0;

    virtual std::vector<vec> record(OutputType) = 0;

    [[nodiscard]] virtual double get_strain_energy() const = 0;
    [[nodiscard]] virtual double get_complementary_energy() const = 0;
    [[nodiscard]] virtual double get_kinetic_energy() const = 0;
    [[nodiscard]] virtual double get_viscous_energy() const = 0;
    [[nodiscard]] virtual double get_nonviscous_energy() const = 0;
    [[nodiscard]] virtual const vec& get_momentum() const = 0;
    [[nodiscard]] virtual double get_momentum_component(DOF) const = 0;

    [[nodiscard]] virtual double get_characteristic_length() const = 0;

    [[nodiscard]] virtual mat compute_shape_function(const mat&, unsigned) const = 0;
};

#endif

//! @}
