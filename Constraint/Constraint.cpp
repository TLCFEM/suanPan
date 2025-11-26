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

#include "Constraint.h"

Constraint::Constraint(const unsigned T, const unsigned AT, uvec&& N, uvec&& D, const unsigned S)
    : ConditionalModifier(T, AT, std::move(N), std::move(D))
    , lagrangian_size(S) {}

const sp_vec& Constraint::get_resistance() const { return resistance; }

const vec& Constraint::get_auxiliary_resistance() const { return auxiliary_resistance; }

const sp_mat& Constraint::get_auxiliary_stiffness() const { return auxiliary_stiffness; }

const sp_mat& Constraint::get_stiffness() const { return stiffness; }

const vec& Constraint::get_auxiliary_load() const { return auxiliary_load; }

void Constraint::set_multiplier_size(const unsigned S) {
    lagrangian_size = S;
    stiffness.reset();
    resistance.reset();
}

unsigned Constraint::get_multiplier_size() const { return lagrangian_size; }
