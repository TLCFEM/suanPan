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

#include "GeneralizedAlphaExplicit.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

void GeneralizedAlphaExplicit::update_parameter(const double NT) { DT = NT; }

int GeneralizedAlphaExplicit::process_load_impl(const bool full) {
    const auto D = get_domain();
    auto& W = D->get_factory();

    const sp_d auto current_time = W->get_current_time();
    const sp_d auto trial_time = W->get_trial_time();

    W->update_trial_time(AF * current_time + (1. - AF) * trial_time);

    const auto code = ExplicitIntegrator::process_load_impl(full);

    W->update_trial_time(trial_time);

    return code;
}

int GeneralizedAlphaExplicit::process_constraint_impl(const bool full) {
    const auto D = get_domain();
    auto& W = D->get_factory();

    const sp_d auto current_time = W->get_current_time();
    const sp_d auto trial_time = W->get_trial_time();

    W->update_trial_time(AF * current_time + (1. - AF) * trial_time);

    const auto code = ExplicitIntegrator::process_constraint_impl(full);

    W->update_trial_time(trial_time);

    return code;
}

int GeneralizedAlphaExplicit::correct_trial_status() {
    const auto D = get_domain();
    auto& W = D->get_factory();

    W->update_trial_displacement_by(B * DT * DT * W->get_trial_acceleration());

    return D->update_trial_status();
}

GeneralizedAlphaExplicit::GeneralizedAlphaExplicit(const unsigned T, const double R)
    : ExplicitIntegrator(T)
    , B((R * R - 5. * R + 10) / 6. / (R - 2.) / (R + 1.))
    , AM((2. * R - 1) / (1. + R))
    , AF(AM - .5) {}

bool GeneralizedAlphaExplicit::has_corrector() const { return true; }

void GeneralizedAlphaExplicit::assemble_resistance() {
    const auto D = get_domain();
    auto& W = D->get_factory();

    auto fa = std::async([&] { D->assemble_resistance(); });
    auto fb = std::async([&] { D->assemble_damping_force(); });
    auto fc = std::async([&] { D->assemble_nonviscous_force(); });
    auto fd = std::async([&] { D->assemble_inertial_force(); });

    fa.get();
    fb.get();
    fc.get();
    fd.get();

    W->set_sushi(W->get_trial_resistance() - AF * W->get_incre_resistance() + W->get_trial_damping_force() - AF * W->get_incre_damping_force() + W->get_trial_nonviscous_force() - AF * W->get_incre_nonviscous_force() + W->get_trial_inertial_force() - AM * W->get_incre_inertial_force());
}

vec GeneralizedAlphaExplicit::get_force_residual() { return ExplicitIntegrator::get_force_residual() / (1. - AM); }

vec GeneralizedAlphaExplicit::get_displacement_residual() { return ExplicitIntegrator::get_displacement_residual() / (1. - AM); }

sp_mat GeneralizedAlphaExplicit::get_reference_load() { return ExplicitIntegrator::get_reference_load() / (1. - AM); }

int GeneralizedAlphaExplicit::update_trial_status(bool) {
    const auto D = get_domain();
    auto& W = D->get_factory();

    W->update_incre_displacement(DT * W->get_current_velocity() + (.5 - B) * DT * DT * W->get_current_acceleration());
    W->update_incre_velocity(DT * W->get_current_acceleration());

    return D->update_trial_status();
}

void GeneralizedAlphaExplicit::print() {
    suanpan_info("An explicit integrator using the Generalized-Alpha algorithm.\n");
}
