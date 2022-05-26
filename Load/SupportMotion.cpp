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

#include "SupportMotion.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Load/Amplitude/Amplitude.h>
#include <Solver/Integrator/Integrator.h>
#include <Step/Step.h>

SupportMotion::SupportMotion(const unsigned T, const unsigned ST, const double L, uvec&& N, const unsigned D, const unsigned AT)
    : Load(T, ST, AT, std::forward<uvec>(N), uvec{D}, L) { enable_displacement_control(); }

SupportMotion::SupportMotion(const unsigned T, const unsigned ST, const double L, uvec&& N, uvec&& D, const unsigned AT)
    : Load(T, ST, AT, std::forward<uvec>(N), std::forward<uvec>(D), L) { enable_displacement_control(); }

int SupportMotion::initialize(const shared_ptr<DomainBase>& D) {
    set_end_step(start_step + 1);

    const auto& W = D->get_factory();

    vector<uword> r_dof, e_dof;

    for(const auto I : W->get_reference_dof()) r_dof.emplace_back(I);

    for(const auto I : get_nodal_active_dof(D)) {
        e_dof.emplace_back(I);
        if(!if_contain(r_dof, e_dof.back()).first) r_dof.emplace_back(e_dof.back());
        else suanpan_warning("more than one displacement load are applied on the same DoF.\n");
    }

    encoding = e_dof;

    W->set_reference_dof(uvec(r_dof));
    W->set_reference_size(static_cast<unsigned>(r_dof.size()));

    return Load::initialize(D);
}

int SupportDisplacement::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();

    trial_settlement.zeros(W->get_size());

    trial_settlement(encoding).fill(pattern * magnitude->get_amplitude(W->get_trial_time()));

    return SUANPAN_SUCCESS;
}

int SupportVelocity::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();
    const auto& G = D->get_current_step()->get_integrator();

    trial_settlement.zeros(W->get_size());

    trial_settlement(encoding) = G->from_total_velocity(pattern * magnitude->get_amplitude(W->get_trial_time()), encoding);

    return SUANPAN_SUCCESS;
}

int SupportAcceleration::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();
    const auto& G = D->get_current_step()->get_integrator();

    trial_settlement.zeros(W->get_size());

    trial_settlement(encoding) = G->from_total_acceleration(pattern * magnitude->get_amplitude(W->get_trial_time()), encoding);

    return SUANPAN_SUCCESS;
}
