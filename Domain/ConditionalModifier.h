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
/**
 * @class ConditionalModifier
 * @brief A ConditionalModifier class.
 *
 * The ConditionalModifier class.
 *
 * TODO: update node/element tags from groups need to be done in each iteration
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

#include <Domain/Tag.h>

class DomainBase;
class Amplitude;

class ConditionalModifier : public Tag {
protected:
    const bool initialized = false;

    const bool connected = false;

    unsigned start_step, end_step = static_cast<unsigned>(-1);

    const unsigned amplitude_tag;

    uvec node_encoding; // node/element encoding
    uvec dof_reference; // reference DoF ZERO based
    uvec dof_encoding;  // DoF encoding

    shared_ptr<Amplitude> magnitude;

    /**
     * \brief Generate active DoF vector from assigned nodes.
     * \return vector of active DoFs
     */
    uvec get_nodal_active_dof(const shared_ptr<DomainBase>&);
    /**
     * \brief Generate active DoF vector from all nodes in the model.
     * \return vector of active DoFs
     */
    uvec get_all_nodal_active_dof(const shared_ptr<DomainBase>&);

public:
    ConditionalModifier(unsigned, unsigned, unsigned, uvec&&, uvec&&);

    virtual int initialize(const shared_ptr<DomainBase>&);

    /**
     * \brief  This method provides all necessary pieces of typical constraints/loads
     * required, including additional blocks in original global stiffness, border matrix
     * resistance of multiplier, external loads.
     * \return success flag
     */
    virtual int process(const shared_ptr<DomainBase>&) = 0;
    /**
     * \brief For some algorithms, the global stiffness is formed only once in each substep.
     * After calling solver, the storage may contain factorization. It is not correct to modify it
     * in those algorithms. This method should provide updated constraint/load resistance but must not
     * touch global stiffness.
     * \return success flag
     */
    virtual int process_resistance(const shared_ptr<DomainBase>&);

    [[nodiscard]] const uvec& get_node_encoding() const;
    [[nodiscard]] const uvec& get_dof_encoding() const;

    void set_initialized(bool) const;
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
    void set_connected(bool) const;
    [[nodiscard]] bool is_connected() const;

    [[nodiscard]] bool validate_step(const shared_ptr<DomainBase>&) const;

    // some may manage state
    virtual void update_status(const vec&);
    virtual void commit_status();
    virtual void clear_status();
    virtual void reset_status();
};

#endif

//! @}
