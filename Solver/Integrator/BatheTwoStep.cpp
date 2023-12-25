/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "BatheTwoStep.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

BatheTwoStep::BatheTwoStep(const unsigned T, const double R, const double G)
    : ImplicitIntegrator(T)
    , GM(G)
    , Q1((R + 1) / (2. * GM * (R - 1) + 4))
    , Q2(.5 - GM * Q1)
    , Q0(1. - Q1 - Q2) {}

void BatheTwoStep::assemble_resistance() {
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

    W->set_sushi(W->get_trial_resistance() + W->get_trial_damping_force() + W->get_trial_nonviscous_force() + W->get_trial_inertial_force());
}

void BatheTwoStep::assemble_matrix() {
    const auto D = get_domain();
    auto& W = D->get_factory();

    auto fa = std::async([&] { D->assemble_trial_stiffness(); });
    auto fb = std::async([&] { D->assemble_trial_geometry(); });
    auto fc = std::async([&] { D->assemble_trial_damping(); });
    auto fd = std::async([&] { D->assemble_trial_nonviscous(); });
    auto fe = std::async([&] { D->assemble_trial_mass(); });

    fa.get();
    fb.get();
    fc.get();
    fd.get();
    fe.get();

    if(W->is_nlgeom()) W->get_stiffness() += W->get_geometry();

    W->get_stiffness() += FLAG::TRAP == step_flag ? P3 * W->get_mass() : P9 * W->get_mass();

    const auto damping_coef = FLAG::TRAP == step_flag ? P2 : P8;

    W->get_stiffness() += W->is_nonviscous() ? damping_coef * (W->get_damping() + W->get_nonviscous()) : damping_coef * W->get_damping();
}

void BatheTwoStep::update_incre_time(double T) {
    const auto& W = get_domain()->get_factory();
    update_parameter(T *= 2.);
    W->update_incre_time(T * (FLAG::TRAP == step_flag ? GM : 1. - GM));
}

int BatheTwoStep::update_trial_status() {
    const auto D = get_domain();

    if(auto& W = D->get_factory(); FLAG::TRAP == step_flag) {
        W->update_trial_acceleration(P3 * W->get_incre_displacement() - P4 * W->get_current_velocity() - W->get_current_acceleration());
        W->update_trial_velocity(P2 * W->get_incre_displacement() - W->get_current_velocity());
    }
    else {
        W->update_trial_velocity(P8 * (W->get_trial_displacement() - W->get_pre_displacement()) - Q02 * W->get_pre_velocity() - Q12 * W->get_current_velocity());
        W->update_trial_acceleration(P8 * (W->get_trial_velocity() - W->get_pre_velocity()) - Q02 * W->get_pre_acceleration() - Q12 * W->get_current_acceleration());
    }

    return D->update_trial_status();
}

void BatheTwoStep::commit_status() {
    const auto D = get_domain();
    auto& W = D->get_factory();

    if(FLAG::TRAP == step_flag) {
        step_flag = FLAG::EULER;
        set_time_step_switch(false);
    }
    else {
        step_flag = FLAG::TRAP;
        set_time_step_switch(true);
    }

    W->commit_pre_displacement();
    W->commit_pre_velocity();
    W->commit_pre_acceleration();

    ImplicitIntegrator::commit_status();
}

void BatheTwoStep::clear_status() {
    step_flag = FLAG::TRAP;
    set_time_step_switch(true);

    ImplicitIntegrator::clear_status();
}

vec BatheTwoStep::from_incre_velocity(const vec& incre_velocity, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    return from_total_velocity(W->get_current_velocity()(encoding) + incre_velocity, encoding);
}

vec BatheTwoStep::from_incre_acceleration(const vec& incre_acceleration, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    return from_total_acceleration(W->get_current_acceleration()(encoding) + incre_acceleration, encoding);
}

vec BatheTwoStep::from_total_velocity(const vec& total_velocity, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    if(FLAG::TRAP == step_flag) return W->get_current_displacement()(encoding) + P1 * (W->get_current_velocity()(encoding) + total_velocity);

    return W->get_pre_displacement()(encoding) + P5 * W->get_pre_velocity()(encoding) + P6 * W->get_current_velocity()(encoding) + P7 * total_velocity;
}

vec BatheTwoStep::from_total_acceleration(const vec& total_acceleration, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    vec total_velocity;
    if(FLAG::TRAP == step_flag) total_velocity = W->get_current_velocity()(encoding) + P1 * (W->get_current_acceleration()(encoding) + total_acceleration);
    else total_velocity = W->get_pre_velocity()(encoding) + P5 * W->get_pre_acceleration()(encoding) + P6 * W->get_current_acceleration()(encoding) + P7 * total_acceleration;

    return from_total_velocity(total_velocity, encoding);
}

void BatheTwoStep::update_parameter(const double NT) {
    if(suanpan::approx_equal(P0, NT)) return;

    P0 = NT;

    P1 = .5 * P0 * GM;
    P2 = 1. / P1;
    P3 = P2 * P2;
    P4 = 2. * P2;

    P5 = P0 * Q0;
    P6 = P0 * Q1;
    P7 = P0 * Q2;
    P8 = 1. / P7;
    P9 = P8 * P8;
}

void BatheTwoStep::print() {
    suanpan_info("A BatheTwoStep time integrator.\n");
}
