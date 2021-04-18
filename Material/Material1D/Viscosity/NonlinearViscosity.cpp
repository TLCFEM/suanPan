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

#include "NonlinearViscosity.h"
#include <Recorder/OutputType.h>

NonlinearViscosity::NonlinearViscosity(const unsigned T, const double A, const double L)
	: DataNonlinearViscosity{fabs(A), fabs(L)}
	, Material1D(T, 0.) {}

void NonlinearViscosity::initialize(const shared_ptr<DomainBase>&) {
	trial_damping = current_damping = initial_damping = compute_damping_coefficient(0., 0.) * (1. == alpha ? 1. : b);

	trial_stiffness = current_stiffness = initial_stiffness.zeros(1);

	trial_strain_rate = current_strain_rate.zeros(1);
}

int NonlinearViscosity::update_trial_status(const vec&) {
	suanpan_error("NonlinearViscosity receives strain only from the associated element, check the model.\n");
	return SUANPAN_FAIL;
}

int NonlinearViscosity::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
	incre_strain = (trial_strain = t_strain) - current_strain;
	incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

	if(norm(incre_strain) + norm(incre_strain_rate) <= 1E-14) return SUANPAN_SUCCESS;

	const auto &u = trial_strain(0), &v = trial_strain_rate(0);

	const auto abs_v = fabs(v);

	const auto eta = compute_damping_coefficient(u, v);

	double term_a, term_b;
	if(0. == alpha) {
		term_a = 0.;
		term_b = 1.;
	}
	else if(1. == alpha) {
		term_a = 1.;
		term_b = v;
	}
	else if(1. < alpha || limit < abs_v) {
		term_a = alpha * pow(abs_v, alpha - 1.);
		term_b = term_a * v / alpha;
	}
	else {
		term_a = 3. * a * v * v + b;
		term_b = (a * v * v + b) * v;
	}

	trial_stress = eta * term_b;
	trial_stiffness = compute_du(u, v) * term_b;
	trial_damping = eta * term_a + compute_dv(u, v) * term_b;

	return SUANPAN_SUCCESS;
}

int NonlinearViscosity::clear_status() {
	current_strain.zeros();
	current_strain_rate.zeros();
	current_stress.zeros();
	current_damping = initial_damping;
	current_stiffness = initial_stiffness;
	return reset_status();
}

int NonlinearViscosity::commit_status() {
	current_strain = trial_strain;
	current_strain_rate = trial_strain_rate;
	current_stress = trial_stress;
	current_damping = trial_damping;
	current_stiffness = trial_stiffness;
	return SUANPAN_SUCCESS;
}

int NonlinearViscosity::reset_status() {
	trial_strain = current_strain;
	trial_strain_rate = current_strain_rate;
	trial_stress = current_stress;
	trial_damping = current_damping;
	trial_stiffness = current_stiffness;
	return SUANPAN_SUCCESS;
}

vector<vec> NonlinearViscosity::record(const OutputType P) {
	vector<vec> data;

	if(OutputType::S == P) data.emplace_back(current_stress);
	else if(OutputType::E == P || OutputType::ED == P) data.emplace_back(current_strain);
	else if(OutputType::V == P || OutputType::VD == P) data.emplace_back(current_strain_rate);

	return data;
}

void NonlinearViscosity::print() {
	suanpan_info("A 1D vicosous damping material %u.\n", get_tag());
	suanpan_info("current strain: %.3E\tcurrent stress: %.3E\n", current_strain(0), current_stress(0));
}
