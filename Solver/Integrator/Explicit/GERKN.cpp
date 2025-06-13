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

#include "GERKN.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

void GERKN::update_parameter(const double NT) { DT = NT; }

int GERKN::process_load_impl(const bool full) {
    auto& W = get_domain()->get_factory();

    const sp_d auto trial_time = W->get_trial_time();

    if(FLAG::FIRST != step_flag) W->update_incre_time((C2 - C1) * DT);

    const auto code = ExplicitIntegrator::process_load_impl(full);

    W->update_trial_time(trial_time);

    return code;
}

int GERKN::process_constraint_impl(const bool full) {
    auto& W = get_domain()->get_factory();

    const sp_d auto trial_time = W->get_trial_time();

    if(FLAG::FIRST != step_flag) W->update_incre_time((C2 - C1) * DT);

    const auto code = ExplicitIntegrator::process_constraint_impl(full);

    W->update_trial_time(trial_time);

    return code;
}

bool GERKN::has_corrector() const { return FLAG::FIRST != step_flag; }

int GERKN::correct_trial_status() {
    // it is guaranteed that it is the second stage when this method is called
    // if(FLAG::FIRST == step_flag) return SUANPAN_SUCCESS;

    const auto D = get_domain();
    auto& W = D->get_factory();

    W->update_trial_displacement(W->get_pre_displacement() + DT * W->get_pre_velocity() + DT * DT * UB0 * W->get_pre_acceleration() + DT * DT * UB1 * W->get_current_acceleration() + DT * DT * UB2 * W->get_trial_acceleration());
    W->update_trial_velocity(W->get_pre_velocity() + DT * VB0 * W->get_pre_acceleration() + DT * VB1 * W->get_current_acceleration() + DT * VB2 * W->get_trial_acceleration());
    W->update_trial_acceleration(AB0 * W->get_pre_acceleration() + AB1 * W->get_current_acceleration() + AB2 * W->get_trial_acceleration());

    return D->update_trial_status();
}

void GERKN::update_incre_time(double T) {
    const auto& W = get_domain()->get_factory();
    update_parameter(T *= 2.);
    W->update_incre_time(T * (FLAG::FIRST == step_flag ? C1 : 1. - C1));
}

int GERKN::update_trial_status(bool) {
    const auto D = get_domain();

    if(auto& W = D->get_factory(); FLAG::FIRST == step_flag) {
        W->update_incre_velocity(DT * VA10 * W->get_current_acceleration());
        W->update_incre_displacement(DT * C1 * W->get_current_velocity() + DT * DT * UA10 * W->get_current_acceleration());
    }
    else {
        W->update_incre_velocity(DT * (VA20 - VA10) * W->get_pre_acceleration() + DT * VA21 * W->get_current_acceleration());
        W->update_incre_displacement(DT * (C2 - C1) * W->get_pre_velocity() + DT * DT * (UA20 - UA10) * W->get_pre_acceleration() + DT * DT * UA21 * W->get_current_acceleration());
    }

    return D->update_trial_status();
}

void GERKN::commit_status() {
    auto& W = get_domain()->get_factory();

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

void GERKN::clear_status() {
    step_flag = FLAG::FIRST;
    set_time_step_switch(true);

    ExplicitIntegrator::clear_status();
}

vec GERKN::from_incre_acceleration(const vec& incre_acceleration, const uvec& encoding) { return from_total_acceleration(get_domain()->get_factory()->get_current_acceleration()(encoding) + incre_acceleration, encoding); }

vec GERKN::from_total_acceleration(const vec& total_acceleration, const uvec& encoding) {
    if(FLAG::FIRST == step_flag) return total_acceleration;

    auto& W = get_domain()->get_factory();

    return total_acceleration / AB2 - AB0 / AB2 * W->get_pre_acceleration()(encoding) - AB1 / AB2 * W->get_current_acceleration()(encoding);
}

WAT2::WAT2(const unsigned T, double R)
    : GERKN(T) {
    R = std::min(1., std::max(0., R));
    R = .1 * R + .3; // [0.3,0.4]

    const auto E = 2. * R - 1.;
    const auto D = 3. * R * R - 3. * R + 1.;
    const auto F = R * D;
    const auto E2 = E * E;

    C1 = R;
    C2 = (3. * C1 - 2.) / 3. / E;

    VA10 = C1;
    VA21 = 2. / 9. * D / C1 / E2;
    VA20 = C2 - VA21;

    UA10 = C1 / 6.;
    UA20 = (3. * C1 - 1.) / 18. / C1;
    UA21 = D / 18. / C1 / E2;

    VB0 = 0.;
    VB1 = .25 / D;
    VB2 = 1. - VB1;

    UB0 = (4. * C1 - 1.) / 12. / C1;
    UB1 = (6. * C1 * C1 - 5. * C1 + 2.) / 24. / F;
    UB2 = .5 - UB0 - UB1;

    AB0 = (1. - 3. * C1) / C1;
    AB1 = (-3. * C1 * C1 + 4. * C1 - 1.) / F;
    AB2 = 1. - AB0 - AB1;
}

void WAT2::print() { suanpan_info("A generalized explicit RKN time integrator using WAT2 scheme with C1={}. doi:10.1002/nme.7658\n", C1); }
