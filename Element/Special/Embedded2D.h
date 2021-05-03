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
 * @class Embedded2D
 * @brief A Embedded2D class.
 *
 * The Embedded2D class.
 *
 * @author tlc
 * @date 12/08/2020
 * @version 0.1.0
 * @file Embedded2D.h
 * @addtogroup Constraint
 * @{
 */

#ifndef EMBEDDED2D_H
#define EMBEDDED2D_H

#include <Element/Element.h>

class Embedded2D final : public Element {
	static constexpr unsigned e_dof = 2;
	static constexpr unsigned max_iteration = 20;

	const unsigned host_tag;
	const unsigned host_size = 0;
	const double alpha;
	const rowvec iso_n;

	uvec idx_x, idx_y;

	shared_ptr<Element> host_element;
public:
	Embedded2D(unsigned, unsigned, unsigned, double = 1E6);

	void initialize(const shared_ptr<DomainBase>&) override;

	int update_status() override;

	int clear_status() override;
	int commit_status() override;
	int reset_status() override;
};

#endif

//! @}
