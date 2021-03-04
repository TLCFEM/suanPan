////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2021 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "GeneralizedAlpha.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

GeneralizedAlpha::GeneralizedAlpha(const unsigned T, const double R)
	: Integrator(T)
	, alpha_f(R / (R + 1.))
	, alpha_m((2. * R - 1.) / (R + 1.))
	, gamma(.5 - (R - 1.) / (R + 1.))
	, beta(pow(R + 1., -2.)) {
	F1 = alpha_f;
	F2 = 1. - F1;
	F3 = alpha_m;
	F4 = 1. - F3;
	F9 = -.5 / beta;
}

GeneralizedAlpha::GeneralizedAlpha(const unsigned T, const double AF, const double AM)
	: Integrator(T)
	, alpha_f(std::min(.5, std::max(AF, .0)))
	, alpha_m(std::min(alpha_f, std::max(AM, -1.)))
	, gamma(.5 - alpha_m + alpha_f)
	, beta(.25 * (gamma + .5) * (gamma + .5)) {
	if(alpha_m != AM || alpha_f != AF) suanpan_error("GeneralizedAlpha() parameters are not acceptable hence automatically adjusted.\n"); // NOLINT(clang-diagnostic-float-equal)

	F1 = alpha_f;
	F2 = 1. - F1;
	F3 = alpha_m;
	F4 = 1. - F3;
	F9 = -.5 / beta;
}

vec GeneralizedAlpha::get_auxiliary_residual() {
	const auto& W = get_domain().lock()->get_factory();

	return W->get_auxiliary_load() - F1 * W->get_current_auxiliary_resistance() - F2 * W->get_trial_auxiliary_resistance();
}

void GeneralizedAlpha::assemble_resistance() {
	const auto& D = get_domain().lock();
	auto& W = D->get_factory();

	D->assemble_resistance();
	D->assemble_inertial_force();
	D->assemble_damping_force();

	W->set_sushi(F1 * (W->get_current_resistance() + W->get_current_damping_force()) + F2 * (W->get_trial_resistance() + W->get_trial_damping_force()) + F3 * W->get_current_inertial_force() + F4 * W->get_trial_inertial_force());
}

void GeneralizedAlpha::assemble_matrix() {
	const auto& D = get_domain().lock();
	auto& W = D->get_factory();

	D->assemble_trial_stiffness();
	D->assemble_trial_mass();
	D->assemble_trial_damping();

	auto& t_stiffness = W->get_stiffness();

	if(W->get_nlgeom()) {
		D->assemble_trial_geometry();
		t_stiffness += W->get_geometry();
	}

	t_stiffness *= F2;
	t_stiffness += F5 * W->get_mass() + F6 * W->get_damping();
}

int GeneralizedAlpha::process_constraint() const {
	const auto& D = get_domain().lock();
	auto& W = D->get_factory();

	const auto current_time = W->get_current_time();
	const auto trial_time = W->get_trial_time();

	W->update_trial_time(F1 * current_time + F2 * trial_time);

	const auto code = D->process_constraint();

	W->update_trial_time(trial_time);

	return code;
}

int GeneralizedAlpha::process_load() const {
	const auto& D = get_domain().lock();
	auto& W = D->get_factory();

	const auto current_time = W->get_current_time();
	const auto trial_time = W->get_trial_time();

	W->update_trial_time(F1 * current_time + F2 * trial_time);

	const auto code = D->process_load();

	W->update_trial_time(trial_time);

	return code;
}

int GeneralizedAlpha::update_trial_status() {
	const auto& D = get_domain().lock();
	auto& W = D->get_factory();

	W->update_incre_acceleration(F7 * W->get_incre_displacement() + F8 * W->get_current_velocity() + F9 * W->get_current_acceleration());
	W->update_incre_velocity(F10 * W->get_current_acceleration() + F11 * W->get_incre_acceleration());

	return D->update_trial_status();
}

void GeneralizedAlpha::update_parameter(const double NT) {
	if(F10 == NT) return; // NOLINT(clang-diagnostic-float-equal)

	F10 = NT;
	F11 = F10 * gamma;
	F8 = -1. / beta / F10;
	F7 = -F8 / F10;
	F6 = -gamma * F8 * F2;
	F5 = F4 * F7;
}

/**
 * \brief update acceleration and velocity for zero displacement increment
 */
void GeneralizedAlpha::update_compatibility() const {
	const auto& D = get_domain().lock();
	auto& W = D->get_factory();

	W->update_incre_acceleration(F8 * W->get_current_velocity() + F9 * W->get_current_acceleration());
	W->update_incre_velocity(F10 * W->get_current_acceleration() + F11 * W->get_incre_acceleration());

	auto& trial_dsp = W->get_trial_displacement();
	auto& trial_vel = W->get_trial_velocity();
	auto& trial_acc = W->get_trial_acceleration();

	auto& t_node_pool = D->get_node_pool();

	suanpan_for_each(t_node_pool.cbegin(), t_node_pool.cend(), [&](const shared_ptr<Node>& t_node) { t_node->update_trial_status(trial_dsp, trial_vel, trial_acc); });
}

void GeneralizedAlpha::print() { suanpan_info("A time integrator using the Generalized-Alpha algorithm.\ndoi:10.1115/1.2900803\n"); }
