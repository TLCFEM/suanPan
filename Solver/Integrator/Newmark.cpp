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

#include "Newmark.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

Newmark::Newmark(const unsigned T, const double A, const double B)
    : Integrator(T)
    , beta(A)
    , gamma(B) {}

void Newmark::assemble_resistance() {
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

void Newmark::assemble_matrix() {
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

    W->get_stiffness() += W->get_geometry() + C0 * W->get_mass() + C1 * W->get_damping();
}

int Newmark::update_trial_status() {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    W->update_incre_acceleration(C0 * W->get_incre_displacement() - C2 * W->get_current_velocity() - C4 * W->get_current_acceleration());
    W->update_incre_velocity(C5 * W->get_current_acceleration() + C3 * W->get_incre_acceleration());

    return D->update_trial_status();
}

/**
 * \brief update acceleration and velocity for zero displacement increment
 */
void Newmark::update_compatibility() const {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    W->update_incre_acceleration(-C2 * W->get_current_velocity() - C4 * W->get_current_acceleration());
    W->update_incre_velocity(C5 * W->get_current_acceleration() + C3 * W->get_incre_acceleration());

    auto& trial_dsp = W->get_trial_displacement();
    auto& trial_vel = W->get_trial_velocity();
    auto& trial_acc = W->get_trial_acceleration();

    auto& t_node_pool = D->get_node_pool();

    suanpan_for_each(t_node_pool.cbegin(), t_node_pool.cend(), [&](const shared_ptr<Node>& t_node) { t_node->update_trial_status(trial_dsp, trial_vel, trial_acc); });
}

vec Newmark::from_incre_velocity(const vec& incre_velocity, const uvec& encoding) {
    auto& W = get_domain().lock()->get_factory();

    return incre_velocity / C1 + C5 * W->get_current_velocity()(encoding) + (C3 * C4 - C5) / C1 * W->get_current_acceleration()(encoding) + W->get_current_displacement()(encoding);
}

vec Newmark::from_incre_acceleration(const vec& incre_acceleration, const uvec& encoding) {
    auto& W = get_domain().lock()->get_factory();

    return incre_acceleration / C0 + C5 * W->get_current_velocity()(encoding) + C4 / C0 * W->get_current_acceleration()(encoding) + W->get_current_displacement()(encoding);
}

void Newmark::update_parameter(const double NT) {
    if(suanpan::approx_equal(C5, NT)) return;

    C5 = NT;
    C2 = 1. / beta / C5;
    C0 = C2 / C5;
    C1 = C2 * gamma;
    C4 = .5 / beta;
    C3 = C5 * gamma;
}

void Newmark::print() { suanpan_info("A Newmark solver.\n"); }
