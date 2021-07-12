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
 * @class BilinearCC
 * @brief The BilinearCC class.
 * 
 * algorithm verified at 26 April 2019 by tlc
 * 
 * @author tlc
 * @date 26/04/2019
 * @version 0.1.0
 * @file BilinearCC.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef BILINEARCC_H
#define BILINEARCC_H

#include "NonlinearCamClay.h"

struct DataBilinearCC {
	const double a, a_slope;
};

class BilinearCC final : DataBilinearCC, public NonlinearCamClay {
	[[nodiscard]] double compute_a(double) const override;
	[[nodiscard]] double compute_da(double) const override;
public:
	BilinearCC(unsigned,   // tag
	           double,     // elastic modulus
	           double,     // poisson's ratio
	           double,     // beta
	           double,     // m
	           double,     // pt
	           double,     // a
	           double,     // a_slope
	           double = 0. // density
	);

	unique_ptr<Material> get_copy() override;

	void print() override;
};

#endif

//! @}
