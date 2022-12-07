/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "OALTS.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

OALTS::OALTS(const unsigned T, const double R)
    : ImplicitIntegrator(T)
    , A1((4. * R - 4.) / (3. - R))
    , A2(-1. - A1)
    , B0(2. / (1. + R) / (3. - R))
    , B1(2. - 2. * B0 + .5 * A1)
    , B2(.5 * A1 + B0)
    , B10(B1 / B0)
    , B20(B2 / B0) { set_time_step_switch(false); }

int OALTS::initialize() {
    if(SUANPAN_SUCCESS != ImplicitIntegrator::initialize()) return SUANPAN_FAIL;

    auto& W = get_domain()->get_factory();

    W->commit_pre_displacement();
    W->commit_pre_velocity();
    W->commit_pre_acceleration();

    return SUANPAN_SUCCESS;
}

void OALTS::assemble_resistance() {
    const auto& D = get_domain();
    auto& W = D->get_factory();

    auto fa = std::async([&] { D->assemble_resistance(); });
    auto fb = std::async([&] { D->assemble_damping_force(); });
    auto fc = std::async([&] { D->assemble_inertial_force(); });

    fa.get();
    fb.get();
    fc.get();

    W->set_sushi(W->get_trial_resistance() + W->get_trial_damping_force() + W->get_trial_inertial_force());
}

void OALTS::assemble_matrix() {
    const auto& D = get_domain();
    auto& W = D->get_factory();

    auto fa = std::async([&] { D->assemble_trial_stiffness(); });
    auto fb = std::async([&] { D->assemble_trial_geometry(); });
    auto fc = std::async([&] { D->assemble_trial_damping(); });
    auto fd = std::async([&] { D->assemble_trial_mass(); });

    fa.get();
    fb.get();
    fc.get();
    fd.get();

    W->get_stiffness() += W->get_geometry() + P1 * P1 * W->get_mass() + P1 * W->get_damping();
}

int OALTS::update_trial_status() {
    const auto& D = get_domain();
    auto& W = D->get_factory();

    W->update_trial_velocity(P1 * W->get_trial_displacement() + P2 * W->get_current_displacement() + P3 * W->get_pre_displacement() - B10 * W->get_current_velocity() - B20 * W->get_pre_velocity());
    W->update_trial_acceleration(P1 * W->get_trial_velocity() + P2 * W->get_current_velocity() + P3 * W->get_pre_velocity() - B10 * W->get_current_acceleration() - B20 * W->get_pre_acceleration());

    return D->update_trial_status();
}

void OALTS::update_parameter(const double NT) {
    if(suanpan::approx_equal(DT, NT)) return;

    DT = NT;

    P1 = 1. / B0 / DT;
    P2 = P1 * A1;
    P3 = P1 * A2;
}

void OALTS::commit_status() {
    auto& W = get_domain()->get_factory();

    W->commit_pre_displacement();
    W->commit_pre_velocity();
    W->commit_pre_acceleration();

    ImplicitIntegrator::commit_status();
}

vec OALTS::from_incre_velocity(const vec& incre_velocity, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    return from_total_velocity(W->get_current_velocity()(encoding) + incre_velocity, encoding);
}

vec OALTS::from_incre_acceleration(const vec& incre_acceleration, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    return from_total_acceleration(W->get_current_acceleration()(encoding) + incre_acceleration, encoding);
}

vec OALTS::from_total_velocity(const vec& total_velocity, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    return total_velocity / P1 - P2 / P1 * W->get_current_displacement()(encoding) - P3 / P1 * W->get_pre_displacement()(encoding) + B10 / P1 * W->get_current_velocity()(encoding) + B20 / P1 * W->get_pre_velocity()(encoding);
}

vec OALTS::from_total_acceleration(const vec& total_acceleration, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    return from_total_velocity(total_acceleration / P1 - P2 / P1 * W->get_current_velocity()(encoding) - P3 / P1 * W->get_pre_velocity()(encoding) + B10 / P1 * W->get_current_acceleration()(encoding) + B20 / P1 * W->get_pre_acceleration()(encoding), encoding);
}

void OALTS::print() { suanpan_info("A time integrator using the OALTS algorithm.\ndoi:10.1002/nme.6188\n"); }
