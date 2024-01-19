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

#include "LeeNewmarkIterative.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

vec LeeNewmarkIterative::update_by_mode_zero(const double mass_coef, const double stiffness_coef) const {
    auto& W = get_domain()->get_factory();

    const auto kernel = current_mass->make_copy();
    kernel += stiffness_coef / mass_coef * current_stiffness;
    const auto damping_force = current_mass * W->get_trial_velocity();
    vec tmp;
    kernel->solve(tmp, damping_force);
    return mass_coef * (damping_force - current_mass * tmp);
}

vec LeeNewmarkIterative::update_by_mode_one(double, double, int) const { return {}; }

vec LeeNewmarkIterative::update_by_mode_two(double, double, int, int) const { return {}; }

vec LeeNewmarkIterative::update_by_mode_three(double, double, double) const { return {}; }

vec LeeNewmarkIterative::update_by_mode_four(double, double, int, int, int, int, double) const { return {}; }

void LeeNewmarkIterative::update_damping_force() const {
    auto& W = get_domain()->get_factory();

    vec summation(W->get_size(), fill::zeros);

    const auto i = [](const double x) { return static_cast<int>(x); };

    for(const auto& [t, p, zeta, omega] : damping_mode) {
        const auto mass_coef = 4. * zeta * omega, stiffness_coef = 4. * zeta / omega;
        switch(t) {
        case Type::T0:
            summation += update_by_mode_zero(mass_coef, stiffness_coef);
            break;
        case Type::T1:
            summation += update_by_mode_one(mass_coef, stiffness_coef, i(p.front()));
            break;
        case Type::T2:
            summation += update_by_mode_two(mass_coef, stiffness_coef, i(p(0)), i(p(1)));
            break;
        case Type::T3:
            summation += update_by_mode_three(mass_coef, stiffness_coef, p.front());
            break;
        case Type::T4:
            summation += update_by_mode_four(mass_coef, stiffness_coef, i(p(0)), i(p(1)), i(p(2)), i(p(3)), p(4));
            break;
        }
    }

    W->update_trial_damping_force_by(summation);
    W->update_sushi_by(summation);
}

LeeNewmarkIterative::LeeNewmarkIterative(const unsigned T, std::vector<Mode>&& M, const double A, const double B)
    : Newmark(T, A, B)
    , damping_mode(std::move(M)) {
    for(auto& [t, p, zeta, omega] : damping_mode)
        switch(t) {
        case Type::T0:
            break;
        case Type::T1:
            if(suanpan::approx_equal(p(0), 0.)) t = Type::T0;
            break;
        case Type::T2:
            if(suanpan::approx_equal(p(0) + p(1), 0.)) t = Type::T0;
            break;
        case Type::T3:
            if(suanpan::approx_equal(p(0), 0.) || p(0) < -1.) t = Type::T0;
            break;
        case Type::T4:
            if(p(4) < -1.) t = Type::T0;
            else if(suanpan::approx_equal(p(4), 0.)) t = suanpan::approx_equal(p(0) + p(1), 0.) ? Type::T0 : Type::T2;
            else if(suanpan::approx_equal(p(0) + p(1) + p(2) + p(3), 0.)) {
                t = Type::T3;
                p = vec{p(4)};
            }
            break;
        }
}

int LeeNewmarkIterative::process_constraint() {
    const auto D = get_domain();
    auto& W = D->get_factory();

    auto& t_mass = W->modify_mass();
    auto& t_stiffness = W->modify_stiffness();

    current_mass.swap(t_mass);
    current_stiffness.swap(t_stiffness);
    if(SUANPAN_SUCCESS != Newmark::process_constraint()) return SUANPAN_FAIL;
    current_mass.swap(t_mass);
    current_stiffness.swap(t_stiffness);

    update_damping_force();

    return Newmark::process_constraint();
}

int LeeNewmarkIterative::process_constraint_resistance() {
    update_damping_force();

    return Newmark::process_constraint_resistance();
}

void LeeNewmarkIterative::assemble_matrix() {
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

    current_mass = W->get_mass()->make_copy();
    current_stiffness = W->get_stiffness()->make_copy();

    W->get_stiffness() += C0 * W->get_mass();

    W->get_stiffness() += W->is_nonviscous() ? C1 * (W->get_damping() + W->get_nonviscous()) : C1 * W->get_damping();
}
