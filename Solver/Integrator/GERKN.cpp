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

bool GERKN::has_corrector() const { return true; }

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

int GERKN::correct_trial_status() {
    if(FLAG::FIRST == step_flag) return SUANPAN_SUCCESS;

    const auto D = get_domain();
    auto& W = D->get_factory();

    W->update_trial_displacement(W->get_pre_displacement() + DT * W->get_pre_velocity() + DT * DT * UB0 * W->get_pre_acceleration() + DT * DT * UB1 * W->get_current_acceleration() + DT * DT * UB2 * W->get_trial_acceleration());
    W->update_trial_velocity(W->get_pre_velocity() + DT * VB0 * W->get_pre_acceleration() + DT * VB1 * W->get_current_acceleration() + DT * VB2 * W->get_trial_acceleration());
    W->update_trial_acceleration(AB0 * W->get_pre_acceleration() + AB1 * W->get_current_acceleration() + AB2 * W->get_trial_acceleration());

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

void GERKN::print() { suanpan_info("A generalized explicit RKN time integrator.\n"); }
