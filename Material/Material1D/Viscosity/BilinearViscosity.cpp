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

#include "BilinearViscosity.h"

double BilinearViscosity::compute_du(double, double) const { return 0.; }

double BilinearViscosity::compute_dv(double, const double strain_rate) const {
	if(const auto abs_v = fabs(strain_rate); abs_v > yield_strain) return hardening;

	return damping;
}

double BilinearViscosity::compute_damping_coefficient(double, const double strain_rate) const {
	const auto abs_v = fabs(strain_rate);

	if(abs_v <= yield_strain) return damping * strain_rate;

	const auto feedback = yield_stress + hardening * (abs_v - yield_strain);

	return strain_rate > 0. ? feedback : -feedback;
}

BilinearViscosity::BilinearViscosity(const unsigned T, const double C, const double S, const double H)
	: DataBilinearViscosity{fabs(C), fabs(S), fabs(C) * H}
	, NonlinearViscosity(T, 0., 1.) {}

unique_ptr<Material> BilinearViscosity::get_copy() { return make_unique<BilinearViscosity>(*this); }
