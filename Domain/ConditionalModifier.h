/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
 * @brief Abstract base for node/element-scoped, step-conditional modifiers (constraints and loads).
 *
 * Design intent
 * - Provide a unified interface for objects that conditionally modify the global system
 *   (stiffness, resistance, external loads, settlements) over analysis steps.
 * - Decouple addressing (nodes/elements/DOFs) from the concrete behavior implemented by derived classes.
 * - Standardize life-cycle, amplitude handling, and step gating across loads and constraints.
 *
 * Targeting and DOFs
 * - Targets are given via tag lists: target_node and target_element.
 * - DOF intent is specified by:
 *   - dof_order: exact DOF layout required by the target (used for validation when non-empty).
 *   - dof_component: unordered subset of DOFs to act upon (falls back to dof_order when empty).
 * - Collected global DOF indices are cached in target_node_dof and target_element_dof:
 *   - collect_node_dof(): from addressed nodes.
 *   - collect_element_dof(): from nodes connected to addressed elements.
 *
 * Validation
 * - Derived classes opt-in to validation via validate_node() / validate_element().
 * - validate_node_impl()/validate_element_impl() ensure targets exist, are active, and match dof_order exactly
 *   when dof_order is provided.
 *
 * Amplitude handling
 * - amplitude_tag identifies an Amplitude resource (via ResourceHolder).
 * - initialize() resolves the amplitude; defaults to Ramp(0) if missing/inactive.
 * - The amplitude start time is aligned to the accumulated time of steps before start_step.
 * - get_amplitude() samples by the current trial time from the Factory.
 *
 * Step gating
 * - start_step and end_step bound the active interval [start_step, end_step) for this modifier.
 * - validate_step() checks if the modifier is active in the current step and itself is active.
 * - Derived classes may adjust these (e.g., displacement-control loads limit to one step).
 *
 * Life-cycle
 * - Construction: record tags and DOF intent.
 * - initialize():
 *   1) resolve amplitude and set start time,
 *   2) optional node/element layout validation,
 *   3) collect and cache target DOFs,
 *   4) mark as initialized.
 * - process(): pure virtual; derived classes assemble stiffness/resistance/loads/settlements.
 * - process_resistance(): defaults to process(); override to avoid touching stiffness when required by algorithms.
 * - stage(): pre-commit hook for predictor–corrector style updates (e.g., restitution).
 * - deinitialize()/is_initialized(): manage lifetime and reusability.
 *
 * Connectivity and queries
 * - is_connected(): whether this modifier should be treated as an “element” for bandwidth/RCM; default false.
 * - get_involving_nodes(): union of explicitly targeted nodes and nodes connected by targeted elements.
 * - get_node_dof(): expose collected node DOFs.
 *
 * State hooks for advanced constraints
 * - update_status(), clear_status(), commit_status(), reset_status(): optional state/multiplier management.
 *
 * Usage in derived classes
 * - Loads (e.g., NodalForce, NodalDisplacement/SupportMotion, BodyForce, LineUDL, NodalAcceleration, ReferenceForce)
 *   employ get_amplitude() and collected DOFs to assemble trial_load, trial_settlement, or reference_load.
 * - Constraints (e.g., BC, FixedLength, NodeLine/Facet, Embed, ParticleCollision, Rigid/RestitutionWallPenalty, MPC)
 *   use collected DOFs to assemble resistance/stiffness and optional auxiliary blocks.
 *
 * Group targeting
 * - GroupModifier translates group tags to concrete node/element tags (DomainBase::flatten_group), allowing
 *   group-based derived classes to reuse the same initialization and DOF collection logic.
 *
 * Concurrency and storage
 * - ConditionalModifier itself is storage-agnostic; derived classes acquire Factory mutexes when modifying global
 *   matrices/vectors under full/sparse schemes.
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
#include <set>

class ConditionalModifier : public UniqueTag {
    const unsigned amplitude_tag;

    bool initialized = false;

    ResourceHolder<Amplitude> amplitude;

    // indicate which DOF components of the target nodes/elements shall be found
    // the modifier itself will apply itself to those active components
    const std::vector<Node::DOF> dof_component;

protected:
    unsigned start_step{1u}, end_step{static_cast<unsigned>(-1)};

    // indicate the order of DOF components
    // used to validate of the problem is compatible
    // the target nodes must have the exact order of DoFs
    // it can be empty such that the order check is not performed
    const std::vector<Node::DOF> dof_order;

    uvec target_node, target_element, target_node_dof;

    [[nodiscard]] double get_amplitude(const shared_ptr<DomainBase>&) const;

    /**
     * \brief Return the DoF components, falls back to the DoF order.
     *
     * When DoF order is given, DoF components can be omitted.
     */
    [[nodiscard]] const std::vector<Node::DOF>& get_dof_component() const;

    [[nodiscard]] bool validate_node(const shared_ptr<DomainBase>&) const;
    [[nodiscard]] bool validate_element(const shared_ptr<DomainBase>&) const;

    [[nodiscard]] uvec collect_node_dof(const shared_ptr<DomainBase>&) const;

public:
    ConditionalModifier(
        unsigned,                 // tag
        unsigned,                 // amplitude tag
        std::vector<Node::DOF>&&, // dof order
        std::vector<Node::DOF>&&  // dof component (unordered)
    );

    virtual int initialize(const shared_ptr<DomainBase>&);

    /**
     * \brief Process and update both stiffness and resistance.
     *
     * This method provides all necessary pieces of typical constraints/loads
     * required, including additional blocks in original global stiffness, border matrix
     * resistance of multiplier, external loads.
     */
    virtual int process(const shared_ptr<DomainBase>&) = 0;
    /**
     * \brief Process and update resistance.
     *
     * For some algorithms, the global stiffness is formed only once in each substep.
     * After calling solver, the storage may contain factorization. It is not correct to modify it
     * in those algorithms. This method should provide updated constraint/load resistance but must not
     * touch global stiffness.
     */
    virtual int process_resistance(const shared_ptr<DomainBase>&);

    /**
     * \brief Some algorithms need to manually modify some variables after solving.
     *
     * Typical example is the predictor--corrector type algorithms.
     * This method is called before committing trial status to perform necessary operations.
     */
    virtual void stage(const shared_ptr<DomainBase>&) {}

    /**
     * \brief Return a set of all nodes involved.
     *
     * Some may define the interaction between nodes and elements.
     * The nodes connected by elements are also found and returned.
     */
    [[nodiscard]] std::set<uword> get_involving_nodes(const shared_ptr<DomainBase>&) const;

    [[nodiscard]] const uvec& get_node_dof() const;

    void deinitialize();
    [[nodiscard]] bool is_initialized() const;

    void set_start_step(unsigned);
    [[nodiscard]] unsigned get_start_step() const;

    void set_end_step(unsigned);
    [[nodiscard]] unsigned get_end_step() const;

    /**
     * \brief Indicate if this modifier can be deemed as an element that needs to account for connectivity.
     *
     * Some constraints may modify global stiffness matrix so that it needs to be treated as an element
     * which may affect bandwidth of banded storage.
     * By calling this method, the RCM reordering algorithm will take this constraint into consideration.
     * Make sure it is properly overridden in the derived classes.
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
