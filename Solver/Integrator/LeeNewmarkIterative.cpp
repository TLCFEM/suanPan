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

void LeeNewmarkIterative::update_damping_force() const {
    auto& W = get_domain()->get_factory();

    vec summation(W->get_size(), fill::zeros);

    for(auto I = 0llu; I < mass_coef.n_elem; ++I) {
        const auto kernel = current_mass->make_copy();
        kernel += stiffness_coef(I) / mass_coef(I) * current_stiffness;
        const auto damping_force = current_mass * W->get_trial_velocity();
        vec tmp;
        kernel->solve(tmp, damping_force);
        summation += mass_coef(I) * (damping_force - current_mass * tmp);
    }

    W->update_trial_damping_force_by(summation);
    W->update_sushi_by(summation);
}

LeeNewmarkIterative::LeeNewmarkIterative(const unsigned T, vec&& X, vec&& F, const double A, const double B)
    : Newmark(T, A, B)
    , mass_coef(4. * X % F)
    , stiffness_coef(4. * X / F) {}

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

int LeeNewmarkIterative::process_constraint() {
    if(SUANPAN_SUCCESS != Newmark::process_constraint()) return SUANPAN_FAIL;

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

    return SUANPAN_SUCCESS;
}

int LeeNewmarkIterative::process_constraint_resistance() {
    if(SUANPAN_SUCCESS != Newmark::process_constraint_resistance()) return SUANPAN_FAIL;

    update_damping_force();

    return SUANPAN_SUCCESS;
}
