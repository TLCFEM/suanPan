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

#include "BatheExplicit.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

BatheExplicit::BatheExplicit(const unsigned T, const double R)
    : ExplicitIntegrator(T)
    , P(2. / (2. + std::sqrt(2. + 2. * R)))
    , Q1((.5 - P) / P / (1. - P))
    , Q2(.5 - P * Q1)
    , Q0(.5 - Q1 - Q2) {}

bool BatheExplicit::has_corrector() const { return true; }

void BatheExplicit::update_incre_time(double T) {
    const auto& W = get_domain()->get_factory();
    update_parameter(T *= 2.);
    W->update_incre_time(T * (FLAG::FIRST == step_flag ? P : 1. - P));
}

int BatheExplicit::update_trial_status(bool) {
    const auto D = get_domain();

    if(auto& W = D->get_factory(); FLAG::FIRST == step_flag) {
        W->update_incre_velocity(A0 * W->get_current_acceleration());
        W->update_incre_displacement(A0 * W->get_current_velocity() + A1 * W->get_current_acceleration());
    }
    else {
        W->update_incre_velocity(A3 * W->get_current_acceleration());
        W->update_incre_displacement(A3 * W->get_current_velocity() + A4 * W->get_current_acceleration());
    }

    return D->update_trial_status();
}

int BatheExplicit::correct_trial_status() {
    const auto D = get_domain();

    if(auto& W = D->get_factory(); FLAG::FIRST == step_flag) W->update_incre_velocity(A2 * W->get_incre_acceleration());
    else W->update_incre_velocity(A5 * W->get_pre_acceleration() + A6 * W->get_current_acceleration() + A7 * W->get_trial_acceleration());

    return D->update_trial_status();
}

void BatheExplicit::commit_status() {
    const auto D = get_domain();
    auto& W = D->get_factory();

    if(FLAG::FIRST == step_flag) {
        step_flag = FLAG::SECOND;
        set_time_step_switch(false);
    }
    else {
        step_flag = FLAG::FIRST;
        set_time_step_switch(true);
    }

    W->commit_pre_displacement();
    W->commit_pre_velocity();
    W->commit_pre_acceleration();

    ExplicitIntegrator::commit_status();
}

void BatheExplicit::clear_status() {
    step_flag = FLAG::FIRST;
    set_time_step_switch(true);

    ExplicitIntegrator::clear_status();
}

void BatheExplicit::update_parameter(const double NT) {
    if(suanpan::approx_equal(DT, NT)) return;

    DT = NT;

    A0 = P * DT;
    A2 = .5 * A0;
    A1 = A0 * A2;
    A3 = DT - A0;
    A4 = .5 * A3 * A3;
    A5 = Q0 * A3;
    A6 = (.5 + Q1) * A3;
    A7 = Q2 * A3;
}

void BatheExplicit::print() { suanpan_info("An explicit Bathe time integrator. doi:10.1016/j.compstruc.2013.06.007\n"); }
