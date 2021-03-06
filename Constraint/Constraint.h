/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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
	static const double multiplier;

	unsigned num_size; // size of multiplier

	uvec nodes; /**< node indices */
	uvec dofs;  /**< DoF indices */

	vec trial_lambda = zeros(num_size);
	vec current_lambda = zeros(num_size);

	friend void set_constraint_multiplier(double);
public:
	Constraint(unsigned, unsigned, unsigned, uvec&&, uvec&&, unsigned = 0);
	Constraint(const Constraint&) = delete;            // copy forbidden
	Constraint(Constraint&&) = delete;                 // move forbidden
	Constraint& operator=(const Constraint&) = delete; // assign forbidden
	Constraint& operator=(Constraint&&) = delete;      // assign forbidden

	~Constraint() override;

	void set_multiplier_size(unsigned);
	[[nodiscard]] unsigned get_multiplier_size() const;

	void update_incre_lambda(const vec&);

	// some constraint may manage state
	void commit_status() override;
	void clear_status() override;
	void reset_status() override;
};

void set_constraint_multiplier(double);

#endif

//! @}
