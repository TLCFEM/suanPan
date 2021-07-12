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
 * @class NonlinearCamClay
 * @brief The NonlinearCamClay class.
 * 
 * A modified Cam-Clay model allows user-defined hardening response.
 * 
 * algorithm verified at 26 April 2019 by tlc
 * 
 * Reference:
 *     1. Computational Methods for Plasticity: Theory and Applications
 *        [10.1002/9780470694626](https://doi.org/10.1002/9780470694626)
 *        Chapter 10 Section 10.1
 * 
 * @author tlc
 * @date 26/04/2019
 * @version 1.0.0
 * @file NonlinearCamClay.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NONLINEARCAMCLAY_H
#define NONLINEARCAMCLAY_H

#include <Material/Material3D/Material3D.h>

struct DataNonlinearCamClay {
	const double elastic_modulus; // elastic modulus
	const double poissons_ratio;  // poisson's ratio
	const double square_beta;     // beta squared
	const double m;               // radus ratio
	const double pt;              // tensile yield hydrostatic stress
};

class NonlinearCamClay : DataNonlinearCamClay, public Material3D {
	static constexpr unsigned max_iteration = 20;
	static const double sqrt_three_two;
	static const mat unit_dev_tensor;

	const double shear = elastic_modulus / (2. + 2. * poissons_ratio);
	const double bulk = elastic_modulus / (3. - 6. * poissons_ratio);
	const double six_shear = 3. * elastic_modulus / (1. + poissons_ratio);
	const double square_m = m * m;

	[[nodiscard]] virtual double compute_a(double) const = 0;
	[[nodiscard]] virtual double compute_da(double) const = 0;
public:
	NonlinearCamClay(unsigned,   // tag
	                 double,     // elastic modulus
	                 double,     // poisson's ratio
	                 double,     // beta
	                 double,     // m
	                 double,     // pt
	                 double = 0. // density
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	[[nodiscard]] double get_parameter(ParameterType) const override;

	int update_trial_status(const vec&) override;

	int clear_status() override;
	int commit_status() override;
	int reset_status() override;

	void print() override;
};

#endif

//! @}
