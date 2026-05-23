/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#include "UDNewmark.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

UDNewmark::UDNewmark(const unsigned T, const double A, const double B, cx_vec&& M, cx_vec&& S)
    : Newmark(T, A, B)
    , m(std::move(M))
    , s(std::move(S))
    , aux_para(accu(m / s).real()) {}

int UDNewmark::initialize() {
    auto& W = get_domain()->get_factory();

    current_nonviscous.zeros(W->get_size(), m.n_elem);

    return Newmark::initialize();
}

void UDNewmark::update_parameter(const double DT) {
    Newmark::update_parameter(DT);

    const cx_vec t_para = 2. / DT + s;
    s_para = (t_para - 2. * s) / t_para;
    m_para = m / t_para;
    accu_para = accu(m_para).real();
}

void UDNewmark::commit_status() {
    current_nonviscous *= diagmat(s_para);
    current_nonviscous += target_field() * m_para.t();

    Newmark::commit_status();
}

void UDNewmark::clear_status() {
    current_nonviscous.zeros();
    Newmark::clear_status();
}

void UDNewmark::print() {
    suanpan_info("A Newmark integrator with the universal damping model.\n");
}

vec UDDNewmark::target_field() const {
    auto& W = get_domain()->get_factory();

    return W->get_current_resistance() + W->get_trial_resistance();
}
void UDDNewmark::assemble_resistance() {
    UDNewmark::assemble_resistance();

    auto& W = get_domain()->get_factory();

    const vec trial_nonviscous = real(current_nonviscous * s_para + accu_para * W->get_current_resistance() + (accu_para - aux_para) * W->get_trial_resistance());

    W->update_trial_nonviscous_force_by(trial_nonviscous);

    W->update_sushi_by(trial_nonviscous);
}

void UDDNewmark::assemble_effective_matrix() {
    auto& W = get_domain()->get_factory();

    if(W->is_nlgeom()) W->get_stiffness() += W->get_geometry();

    W->get_stiffness() *= 1. + accu_para - aux_para;

    W->get_stiffness() += C0 * W->get_mass();

    W->get_stiffness() += W->is_nonviscous() ? C1 * (W->get_damping() + W->get_nonviscous()) : C1 * W->get_damping();
}

vec UDANewmark::target_field() const {
    auto& W = get_domain()->get_factory();

    return W->get_current_inertial_force() + W->get_trial_inertial_force();
}
void UDANewmark::assemble_resistance() {
    UDNewmark::assemble_resistance();

    auto& W = get_domain()->get_factory();

    const vec trial_nonviscous = real(current_nonviscous * s_para + accu_para * W->get_current_inertial_force() + (accu_para - aux_para) * W->get_trial_inertial_force());

    W->update_trial_nonviscous_force_by(trial_nonviscous);

    W->update_sushi_by(trial_nonviscous);
}

void UDANewmark::assemble_effective_matrix() {
    auto& W = get_domain()->get_factory();

    if(W->is_nlgeom()) W->get_stiffness() += W->get_geometry();

    W->get_stiffness() += (C0 + C0 * (accu_para - aux_para)) * W->get_mass();

    W->get_stiffness() += W->is_nonviscous() ? C1 * (W->get_damping() + W->get_nonviscous()) : C1 * W->get_damping();
}
