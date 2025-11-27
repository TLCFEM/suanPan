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
 * @class ConditionalModifier
 * @brief A ConditionalModifier class.
 *
 * A parent abstract class for constraints and loads.
 * Both are node based and conditionally applied to the system.
 * Thus, they conditionally modify the system.
 *
 * @author tlc
 * @date 07/03/2021
 * @version 0.1.0
 * @file ConditionalModifier.h
 * @addtogroup ConditionalModifier
 * @{
 */

#ifndef CONDITIONALMODIFIER_H
#define CONDITIONALMODIFIER_H

#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Load/Amplitude/Amplitude.h>
#include <Toolbox/ResourceHolder.h>

class ConditionalModifier : public UniqueTag {
    const unsigned amplitude_tag;

    bool initialized = false;

    ResourceHolder<Amplitude> amplitude;

    // indicate which DOF components of the target nodes/elements shall be found
    // the modifier itself will apply itself to those active components
    const std::vector<Node::DOF> dof_component;

    bool validate_dof(const shared_ptr<DomainBase>&);
    uvec update_active_dof(const shared_ptr<DomainBase>&);

protected:
    unsigned start_step{1u}, end_step{static_cast<unsigned>(-1)};

    // indicate the order of DOF components
    // used to validate of the problem is compatible
    // the target nodes must have the exact order of DoFs
    // it can be empty such that the order check is not performed
    const std::vector<Node::DOF> dof_order;

    uvec target_node, target_dof;

    [[nodiscard]] double get_amplitude(const shared_ptr<DomainBase>&) const;

    const std::vector<Node::DOF>& get_dof_component() const;

public:
    ConditionalModifier(
        unsigned,                 // tag
        unsigned,                 // amplitude tag
        uvec&&,                   // object tag
        std::vector<Node::DOF>&&, // dof order
        std::vector<Node::DOF>&&  // dof component (unordered)
    );

    virtual int initialize(const shared_ptr<DomainBase>&);

    /**
     * \brief  This method provides all necessary pieces of typical constraints/loads
     * required, including additional blocks in original global stiffness, border matrix
     * resistance of multiplier, external loads.
     */
    virtual int process(const shared_ptr<DomainBase>&) = 0;
    /**
     * \brief For some algorithms, the global stiffness is formed only once in each substep.
     * After calling solver, the storage may contain factorization. It is not correct to modify it
     * in those algorithms. This method should provide updated constraint/load resistance but must not
     * touch global stiffness.
     */
    virtual int process_resistance(const shared_ptr<DomainBase>&);

    /**
     * \brief Some algorithms need to manually modify some variables after solving. Typical example is the
     * predictor--corrector type algorithms. This method is called before committing trial status to perform
     * necessary operations.
     */
    virtual void stage(const shared_ptr<DomainBase>&) {}

    [[nodiscard]] const uvec& get_node_encoding() const;
    [[nodiscard]] const uvec& get_dof_encoding() const;

    void deinitialize();
    [[nodiscard]] bool is_initialized() const;

    void set_start_step(unsigned);
    [[nodiscard]] unsigned get_start_step() const;

    void set_end_step(unsigned);
    [[nodiscard]] unsigned get_end_step() const;

    /**
     * \brief Some constraints may modify global stiffness matrix so that it needs to be treated as an element
     * which may affect bandwidth of banded storage. By calling this method, the RCM reordering algorithm will
     * take this constraint into consideration. Make sure it is called in the constructor.
     */
    [[nodiscard]] virtual bool is_connected() const { return false; }

    /**
     * \brief Validate itself against the current active step to see if itself needs to be applied.
     */
    [[nodiscard]] bool validate_step(const shared_ptr<DomainBase>&) const;

    // some may manage state
    virtual void update_status(const vec&) {}
    virtual void clear_status() {}
    virtual void commit_status() {}
    virtual void reset_status() {}
};

class GroupModifier {
    const uvec groups;

protected:
    [[nodiscard]] uvec update_object_tag(const shared_ptr<DomainBase>&) const;

public:
    explicit GroupModifier(uvec&&);
};

std::vector<Node::DOF> parse_dof(std::string_view);

#endif

//! @}
