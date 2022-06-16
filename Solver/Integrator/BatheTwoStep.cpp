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

#include "BatheTwoStep.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

void BatheTwoStep::assemble_resistance() {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    auto fa = std::async([&] { D->assemble_resistance(); });
    auto fb = std::async([&] { D->assemble_damping_force(); });
    auto fc = std::async([&] { D->assemble_inertial_force(); });

    fa.get();
    fb.get();
    fc.get();

    W->set_sushi(W->get_trial_resistance() + W->get_trial_damping_force() + W->get_trial_inertial_force());
}

void BatheTwoStep::assemble_matrix() {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    auto fa = std::async([&] { D->assemble_trial_stiffness(); });
    auto fb = std::async([&] { D->assemble_trial_geometry(); });
    auto fc = std::async([&] { D->assemble_trial_damping(); });
    auto fd = std::async([&] { D->assemble_trial_mass(); });

    fa.get();
    fb.get();
    fc.get();
    fd.get();

    auto& t_stiff = W->get_stiffness();

    t_stiff += W->get_geometry();

    t_stiff += FLAG::TRAP == step_flag ? C6 * W->get_mass() + C3 * W->get_damping() : C5 * W->get_mass() + C2 * W->get_damping();
}

int BatheTwoStep::update_trial_status() {
    const auto& D = get_domain().lock();

    if(auto& W = D->get_factory(); FLAG::TRAP == step_flag) {
        W->update_incre_acceleration(C6 * W->get_incre_displacement() - C4 * W->get_current_velocity() - 2. * W->get_current_acceleration());
        W->update_incre_velocity(C3 * W->get_incre_displacement() - 2. * W->get_current_velocity());
    }
    else {
        W->update_trial_velocity(C2 * W->get_incre_displacement() + C1 * (W->get_pre_displacement() - W->get_current_displacement()));
        W->update_trial_acceleration(C1 * W->get_pre_velocity() - C3 * W->get_current_velocity() + C2 * W->get_trial_velocity());
    }

    return D->update_trial_status();
}

void BatheTwoStep::commit_status() {
    const auto& D = get_domain().lock();
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

    Integrator::commit_status();
}

void BatheTwoStep::clear_status() {
    step_flag = FLAG::TRAP;
    set_time_step_switch(true);

    Integrator::clear_status();
}

/**
 * \brief update acceleration and velocity for zero displacement increment
 */
void BatheTwoStep::update_compatibility() const {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    if(FLAG::TRAP == step_flag) {
        W->update_incre_acceleration(-C4 * W->get_current_velocity() - 2. * W->get_current_acceleration());
        W->update_incre_velocity(-2. * W->get_current_velocity());
    }
    else {
        W->update_trial_velocity(C1 * (W->get_pre_displacement() - W->get_current_displacement()));
        W->update_trial_acceleration(C1 * W->get_pre_velocity() - C3 * W->get_current_velocity() + C2 * W->get_trial_velocity());
    }

    auto& trial_dsp = W->get_trial_displacement();
    auto& trial_vel = W->get_trial_velocity();
    auto& trial_acc = W->get_trial_acceleration();

    suanpan::for_all(D->get_node_pool(), [&](const shared_ptr<Node>& t_node) { t_node->update_trial_status(trial_dsp, trial_vel, trial_acc); });
}

vec BatheTwoStep::from_incre_velocity(const vec& incre_velocity, const uvec& encoding) {
    auto& W = get_domain().lock()->get_factory();

    vec total_displacement = W->get_current_displacement()(encoding);

    if(FLAG::TRAP == step_flag) total_displacement += incre_velocity / C3 + C0 * W->get_current_velocity()(encoding);
    else total_displacement += (incre_velocity + W->get_current_velocity()(encoding)) / C2 + (W->get_current_displacement()(encoding) - W->get_pre_displacement()(encoding)) / 3.;

    return total_displacement;
}

vec BatheTwoStep::from_incre_acceleration(const vec& incre_acceleration, const uvec& encoding) {
    auto& W = get_domain().lock()->get_factory();

    vec total_displacement = W->get_current_displacement()(encoding);

    if(FLAG::TRAP == step_flag) total_displacement += incre_acceleration / C6 + C0 * W->get_current_velocity()(encoding) + 2. / C6 * W->get_current_acceleration()(encoding);
    else total_displacement += (incre_acceleration + W->get_current_acceleration()(encoding)) / C5 + C3 / C5 * W->get_current_velocity()(encoding) - C1 / C5 * W->get_pre_velocity()(encoding) + (W->get_current_displacement()(encoding) - W->get_pre_displacement()(encoding)) / 3.;

    return total_displacement;
}

void BatheTwoStep::update_parameter(const double NT) {
    if(suanpan::approx_equal(C0, NT)) return;

    C0 = NT;
    C1 = .5 / C0;
    C2 = 3. * C1;
    C3 = 4. * C1;
    C4 = 2. * C3;
    C5 = C2 * C2;
    C6 = C4 / C0;
}

void BatheTwoStep::print() { suanpan_info("A BatheTwoStep solver.\n"); }
