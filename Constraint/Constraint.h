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
 * @class Constraint
 * @brief A Constraint class.
 *
 * The Constraint class.
 *
 * @author tlc
 * @date 03/07/2017
 * @version 0.1.0
 * @file Constraint.h
 * @addtogroup Constraint
 * @{
 */

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <Domain/ConditionalModifier.h>

class Constraint : public ConditionalModifier {
protected:
    unsigned lagrangian_size; // size of multiplier

    vec trial_lambda{lagrangian_size, fill::zeros};
    vec current_lambda{lagrangian_size, fill::zeros};

    sp_vec resistance;
    sp_mat stiffness;

    vec auxiliary_resistance;
    vec auxiliary_load;
    sp_mat auxiliary_stiffness;

public:
    Constraint(
        unsigned,                 // tag
        unsigned,                 // amplitude tag
        uvec&&,                   // node tags
        std::set<Node::DOF>&&,    // dof component (unordered)
        std::vector<Node::DOF>&&, // dof order
        unsigned                  // size of multiplier
    );

    const sp_vec& get_resistance() const;
    const sp_mat& get_stiffness() const;

    const vec& get_auxiliary_resistance() const;
    const vec& get_auxiliary_load() const;
    const sp_mat& get_auxiliary_stiffness() const;

    void set_multiplier_size(unsigned);
    [[nodiscard]] unsigned get_multiplier_size() const;
};

#endif

//! @}
