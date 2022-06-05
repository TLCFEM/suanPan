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

#include "GSSSS.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

GSSSS::GSSSS(const unsigned T)
    : Integrator(T) {}

void GSSSS::assemble_resistance() {
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

void GSSSS::assemble_matrix() {
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

    W->get_stiffness() += W->get_geometry() + XCVD * W->get_damping() + XCAD * W->get_mass();
}

int GSSSS::process_load() {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    const auto current_time = W->get_current_time();
    const auto trial_time = W->get_trial_time();

    W->update_trial_time((1. - W1) * current_time + W1 * trial_time);

    const auto code = Integrator::process_load();

    W->update_trial_time(trial_time);

    return code;
}

int GSSSS::process_constraint() {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    const auto current_time = W->get_current_time();
    const auto trial_time = W->get_trial_time();

    W->update_trial_time((1. - W1) * current_time + W1 * trial_time);

    const auto code = Integrator::process_constraint();

    W->update_trial_time(trial_time);

    return code;
}

int GSSSS::process_load_resistance() {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    const auto current_time = W->get_current_time();
    const auto trial_time = W->get_trial_time();

    W->update_trial_time((1. - W1) * current_time + W1 * trial_time);

    const auto code = Integrator::process_load_resistance();

    W->update_trial_time(trial_time);

    return code;
}

int GSSSS::process_constraint_resistance() {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    const auto current_time = W->get_current_time();
    const auto trial_time = W->get_trial_time();

    W->update_trial_time((1. - W1) * current_time + W1 * trial_time);

    const auto code = Integrator::process_constraint_resistance();

    W->update_trial_time(trial_time);

    return code;
}

int GSSSS::update_trial_status() {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    W->update_trial_velocity(W->get_trial_velocity() + XCVD * W->get_ninja());
    W->update_trial_acceleration(W->get_trial_acceleration() + XCAD * W->get_ninja());

    return D->update_trial_status();
}

void GSSSS::stage_status() {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    W->update_incre_acceleration(W->get_incre_acceleration() / W1G6);
    W->update_incre_velocity(L4 * DT * W->get_current_acceleration() + L5 * DT * W->get_incre_acceleration());
    W->update_incre_displacement(L1 * DT * W->get_current_velocity() + L2 * DT * DT * W->get_current_acceleration() + L3 * DT * DT * W->get_incre_acceleration());

    // since iterative result does not equal to committed result
    // need to sync with elements and nodes
    [[maybe_unused]] const auto code = D->update_trial_status();

    Integrator::stage_status();
}

void GSSSS::update_parameter(const double NT) {
    if(suanpan::approx_equal(DT, NT)) return;

    DT = NT;

    XPV3 = W1G4 * DT - W2G5 * L2 / L3 * DT;
    XPA2 = -W1G6 * L1 / L3 / DT;
    XCVD = W2G5 / DT / W3G3;
    XCAD = W1G6 / DT / DT / W3G3;
}

void GSSSS::update_compatibility() const {
    const auto& D = get_domain().lock();
    auto& W = D->get_factory();

    auto fb = std::async([&] { W->update_trial_velocity(XPV2 * W->get_current_velocity() + XPV3 * W->get_current_acceleration()); });
    auto fc = std::async([&] { W->update_trial_acceleration(XPA2 * W->get_current_velocity() + XPA3 * W->get_current_acceleration()); });

    fb.get();
    fc.get();

    auto& trial_dsp = W->get_trial_displacement();
    auto& trial_vel = W->get_trial_velocity();
    auto& trial_acc = W->get_trial_acceleration();

    auto& t_node_pool = D->get_node_pool();

    suanpan_for_each(t_node_pool.cbegin(), t_node_pool.cend(), [&](const shared_ptr<Node>& t_node) { t_node->update_trial_status(trial_dsp, trial_vel, trial_acc); });
}

vec GSSSS::from_incre_velocity(const vec& incre_velocity, const uvec& encoding) {
    auto& W = get_domain().lock()->get_factory();

    return W->get_current_displacement()(encoding) + L1 * DT * W->get_current_velocity()(encoding) + (L2 - L3 * L4 / L5) * DT * DT * W->get_current_acceleration()(encoding) + L3 / L5 * DT * incre_velocity;
}

vec GSSSS::from_incre_acceleration(const vec& incre_acceleration, const uvec& encoding) {
    auto& W = get_domain().lock()->get_factory();

    return W->get_current_displacement()(encoding) + L1 * DT * W->get_current_velocity()(encoding) + L2 * DT * DT * W->get_current_acceleration()(encoding) + L3 * DT * DT * incre_acceleration;
}

void GSSSS::print() { suanpan_info("A time integrator using the GSSSS algorithm.\n"); }

template<> void GSSSS::generate_constants<GSSSSU0>(const double R3, const double R1, const double R2) {
    L1 = 1.;
    L2 = .5;
    L3 = 1. / (1. + R1) / (1. + R2);
    L4 = 1.;
    L5 = .5 * (3. + R1 + R2 - R1 * R2) * L3;

    W1G1 = 1. / (1. + R3);
    W2G2 = .5 * W1G1;
    W3G3 = L3 * W1G1;
    W1G4 = W1G1;
    W2G5 = L5 * W1G1;
    W1G6 = (2. + R1 + R2 + R3 - R1 * R2 * R3) * W3G3;

    W1 = W1G1;

    XPV2 = 1. - W2G5 * L1 / L3;
    XPA3 = 1. - W1G6 * L2 / L3;
}

GSSSSU0::GSSSSU0(const unsigned T, vec&& R)
    : GSSSS(T) {
    R = sort(R.clamp(0., 1.));
    generate_constants<GSSSSU0>(R(0), R(1), R(2));
}

template<> void GSSSS::generate_constants<GSSSSV0>(const double R3, const double R1, const double R2) {
    L1 = 1.;
    L2 = .5;
    L3 = .5 / (1. + R3);
    L4 = 1.;
    L5 = 2. * L3;

    W2G2 = 1. / (1. + R1) / (1. + R2);
    W3G3 = W2G2 / (1. + R3);
    W1G1 = .5 * (3. + R1 + R2 - R1 * R2) * W2G2;
    W1G4 = W1G1;
    W2G5 = 2. * W3G3;
    W1G6 = (2. + R1 + R2 + R3 - R1 * R2 * R3) * W3G3;

    const auto T0 = 9. - 11. * R1 - 11 * R2 + 19. * R1 * R2;
    const auto T1 = -30. * (3. - 4. * R1 - 4. * R2 + 6. * R1 * R2);
    const auto T2 = 7.5 * (25. - 37. * R1 - 37. * R2 + 53. * R1 * R2);
    const auto T3 = -35. * (3. - 5. * R1 - 5. * R2 + 7. * R1 * R2);

    W1 = (T0 / 2. + T1 / 3. + T2 / 4. + T3 / 5.) / (T0 + T1 / 2. + T2 / 3. + T3 / 4.);

    XPV2 = 1. - W2G5 * L1 / L3;
    XPA3 = 1. - W1G6 * L2 / L3;
}

GSSSSV0::GSSSSV0(const unsigned T, vec&& R)
    : GSSSS(T) {
    R = sort(R.clamp(0., 1.));
    generate_constants<GSSSSV0>(R(0), R(1), R(2));
}
