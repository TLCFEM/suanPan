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
 * @class FixedLength
 * @brief A FixedLength class.
 *
 * @author tlc
 * @date 07/03/2021
 * @version 0.1.0
 * @file FixedLength.h
 * @addtogroup Constraint
 * @{
 */

#ifndef FIXEDLENGTH_H
#define FIXEDLENGTH_H

#include <Constraint/Constraint.h>

class Node;

class FixedLength final : public Constraint {
	weak_ptr<Node> node_i, node_j;
public:
	FixedLength(unsigned, unsigned, unsigned, unsigned, uvec&&);

	int initialize(const shared_ptr<DomainBase>&) override;

	int process(const shared_ptr<DomainBase>&) override;

	void update_status(const vec&) override;
	void commit_status() override;
	void clear_status() override;
	void reset_status() override;
};

#endif

//! @}
