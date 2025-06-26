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

#include "GSSE.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

void GSSE::update_parameter(const double NT) { DT = NT; }

int GSSE::process_load_impl(const bool full) {
    auto& W = get_domain()->get_factory();

    const sp_d auto trial_time = W->get_trial_time();

    W->update_incre_time(C * DT);

    const auto code = ExplicitIntegrator::process_load_impl(full);

    W->update_trial_time(trial_time);

    return code;
}

int GSSE::process_constraint_impl(const bool full) {
    auto& W = get_domain()->get_factory();

    const sp_d auto trial_time = W->get_trial_time();

    W->update_incre_time(C * DT);

    const auto code = ExplicitIntegrator::process_constraint_impl(full);

    W->update_trial_time(trial_time);

    return code;
}

bool GSSE::has_corrector() const { return true; }

int GSSE::correct_trial_status() {
    const auto D = get_domain();
    auto& W = D->get_factory();

    W->update_incre_displacement(DT * W->get_current_velocity() + DT * DT * (.5 - UB) * W->get_current_acceleration() + DT * DT * UB * W->get_trial_acceleration());
    W->update_incre_velocity(DT * (1. - VB) * W->get_current_acceleration() + DT * VB * W->get_trial_acceleration());
    W->update_incre_acceleration(AB * (W->get_trial_acceleration() - W->get_current_acceleration()));

    return D->update_trial_status();
}

GSSE::GSSE(const unsigned T, const double RS, const double RB)
    : ExplicitIntegrator(T)
    , C((2. + RS - RB * RS) / (RB * RS + RB + RS + 1.))
    , UA0(.5 * C * C)
    , VA0(C)
    , UB((RB * (RB - 1) * std::pow(RS, 3) + .5 * (RB - 1.) * (RB * RB + 6. * RB - 1.) * RS * RS + (1. - 5. * RB) * RS + 1. - RB) / (-1. - RB) / (1. + RS) * std::pow(RB * RS - RS - 2., -2))
    , VB(.5 / C)
    , AB(1. / C) {}

GSSE::GSSE(const unsigned T, const double R)
    : GSSE(T, R, R) {}

int GSSE::update_trial_status(bool) {
    const auto D = get_domain();
    auto& W = D->get_factory();

    W->update_incre_velocity(DT * VA0 * W->get_current_acceleration());
    W->update_incre_displacement(DT * C * W->get_current_velocity() + DT * DT * UA0 * W->get_current_acceleration());

    return D->update_trial_status();
}

vec GSSE::from_incre_acceleration(const vec& incre_acceleration, const uvec& encoding) { return 1. / C / AB * incre_acceleration + get_domain()->get_factory()->get_current_acceleration()(encoding); }

vec GSSE::from_total_acceleration(const vec& total_acceleration, const uvec& encoding) { return from_incre_acceleration(total_acceleration - get_domain()->get_factory()->get_current_acceleration()(encoding), encoding); }

void GSSE::print() { suanpan_info("A generalized explicit time integrator. doi:10.1002/nme.6574\n"); }

ICL::ICL(const unsigned T, const double R)
    : GSSE(T, .5 * (1. / R - 1.), R) {}
