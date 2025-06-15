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

#include "GSSSS.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

void GSSSS::update_parameter(const double NT) {
    if(suanpan::approx_equal(DT, NT)) return;

    DT = NT;

    const auto L3T = 1. / L3 / DT;

    C0 = L3T / DT;
    C1 = -L1 * L3T;
    C2 = -L2 / L3;
    C3 = L4 * DT;
    C4 = L5 * DT;

    XD = L3 / W3G3;
    XV = W2G5 / W3G3 / DT;
    XA = W1G6 / W3G3 / DT / DT;
}

int GSSSS::process_load_impl(const bool full) {
    const auto D = get_domain();
    auto& W = D->get_factory();

    const sp_d auto current_time = W->get_current_time();
    const sp_d auto trial_time = W->get_trial_time();

    W->update_trial_time((1. - W1) * current_time + W1 * trial_time);

    const auto code = ImplicitIntegrator::process_load_impl(full);

    W->update_trial_time(trial_time);

    return code;
}

int GSSSS::process_constraint_impl(const bool full) {
    const auto D = get_domain();
    auto& W = D->get_factory();

    const sp_d auto current_time = W->get_current_time();
    const sp_d auto trial_time = W->get_trial_time();

    W->update_trial_time((1. - W1) * current_time + W1 * trial_time);

    const auto code = ImplicitIntegrator::process_constraint_impl(full);

    W->update_trial_time(trial_time);

    return code;
}

GSSSS::GSSSS(const unsigned T)
    : ImplicitIntegrator(T)
    , L1(1.)
    , L2(.5)
    , L4(1.) {}

void GSSSS::assemble_resistance() {
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

    W->set_sushi(W->get_current_resistance() + W3G3 / L3 * W->get_incre_resistance() + W->get_current_damping_force() + W2G5 / L5 * W->get_incre_damping_force() + W->get_current_nonviscous_force() + W2G5 / L5 * W->get_incre_nonviscous_force() + W->get_current_inertial_force() + W1G6 * W->get_incre_inertial_force());
}

void GSSSS::assemble_matrix() {
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

    W->get_stiffness() += XA * W->get_mass();

    W->get_stiffness() += W->is_nonviscous() ? XV * (W->get_damping() + W->get_nonviscous()) : XV * W->get_damping();
}

vec GSSSS::get_force_residual() { return XD * ImplicitIntegrator::get_force_residual(); }

vec GSSSS::get_displacement_residual() { return XD * ImplicitIntegrator::get_displacement_residual(); }

sp_mat GSSSS::get_reference_load() { return XD * ImplicitIntegrator::get_reference_load(); }

int GSSSS::update_trial_status(bool) {
    const auto D = get_domain();
    auto& W = D->get_factory();

    W->update_incre_acceleration(C0 * W->get_incre_displacement() + C1 * W->get_current_velocity() + C2 * W->get_current_acceleration());
    W->update_incre_velocity(C3 * W->get_current_acceleration() + C4 * W->get_incre_acceleration());

    return D->update_trial_status();
}

vec GSSSS::from_incre_velocity(const vec& incre_velocity, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    return from_incre_acceleration(incre_velocity / C4 - W1 * C3 / C4 * W->get_current_acceleration()(encoding), encoding);
}

vec GSSSS::from_incre_acceleration(const vec& incre_acceleration, const uvec& encoding) {
    auto& W = get_domain()->get_factory();

    return incre_acceleration / W1 / C0 - C1 / C0 * W->get_current_velocity()(encoding) - C2 / C0 * W->get_current_acceleration()(encoding) + W->get_current_displacement()(encoding);
}

void GSSSS::print() {
    suanpan_info("A time integrator using the GSSSS algorithm.\n");
}

template<> void GSSSS::generate_constants<GSSSSU0>(const double R3, const double R1, const double R2) {
    L3 = 1. / (1. + R1) / (1. + R2);
    L5 = .5 * (3. + R1 + R2 - R1 * R2) * L3;

    W1 = 1. / (1. + R3);

    W3G3 = L3 * W1;
    W2G5 = L5 * W1;
    W1G6 = (2. + R1 + R2 + R3 - R1 * R2 * R3) * W3G3;
}

GSSSSU0::GSSSSU0(const unsigned T, vec&& R)
    : GSSSS(T) {
    R = sort(R.clamp(0., 1.));
    generate_constants<GSSSSU0>(R(0), R(1), R(2));
}

template<> void GSSSS::generate_constants<GSSSSV0>(const double R3, const double R1, const double R2) {
    L3 = .5 / (1. + R3);
    L5 = 2. * L3;

    const auto T0 = 9. - 11. * R1 - 11. * R2 + 19. * R1 * R2;
    const auto T1 = -30. * (3. - 4. * R1 - 4. * R2 + 6. * R1 * R2);
    const auto T2 = 7.5 * (25. - 37. * R1 - 37. * R2 + 53. * R1 * R2);
    const auto T3 = -35. * (3. - 5. * R1 - 5. * R2 + 7. * R1 * R2);

    W1 = (T0 / 2. + T1 / 3. + T2 / 4. + T3 / 5.) / (T0 + T1 / 2. + T2 / 3. + T3 / 4.);

    W3G3 = 1. / (1. + R1) / (1. + R2) / (1. + R3);
    W2G5 = 2. * W3G3;
    W1G6 = (2. + R1 + R2 + R3 - R1 * R2 * R3) * W3G3;
}

GSSSSV0::GSSSSV0(const unsigned T, vec&& R)
    : GSSSS(T) {
    R = sort(R.clamp(0., 1.));
    generate_constants<GSSSSV0>(R(0), R(1), R(2));
}

template<> void GSSSS::generate_constants<GSSSSOptimal>(const double R, double, double) {
    L3 = .5 / (1. + R);
    L5 = 2. * L3;

    W1 = L5;

    W3G3 = L3 * L5;
    W2G5 = L5 * L5;
    W1G6 = (1. + R) * (3. - R) * W3G3;
}

GSSSSOptimal::GSSSSOptimal(const unsigned T, double R)
    : GSSSS(T) {
    R = std::min(1., std::max(0., R));
    generate_constants<GSSSSOptimal>(R, R, 1.);
}
