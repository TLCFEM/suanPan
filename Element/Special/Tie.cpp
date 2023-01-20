/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include "Tie.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

Tie::Tie(const unsigned T, uvec&& N, uvec&& D, vec&& W, const double L, const double P)
    : Element(T, static_cast<unsigned>(N.size()), static_cast<unsigned>(D.max()), std::forward<uvec>(N), {})
    , dof_pool(std::forward<uvec>(D))
    , weight_pool(std::forward<vec>(W))
    , pseudo_load(L)
    , penalty(P) {}

int Tie::initialize(const shared_ptr<DomainBase>&) {
    initial_stiffness.zeros(get_total_number(), get_total_number());

    const auto n_dof = get_dof_number();

    t_span.zeros(dof_pool.n_elem);
    for(uword I = 0, J = 0; I < dof_pool.n_elem; ++I, J += n_dof) t_span(I) = J + dof_pool(I) - 1;

    initial_stiffness(t_span, t_span) = penalty * weight_pool * weight_pool.t();
    trial_stiffness = current_stiffness = initial_stiffness;

    trial_resistance.zeros(get_total_number());
    trial_resistance(t_span) -= pseudo_load * penalty * weight_pool;

    current_resistance = trial_resistance;

    return SUANPAN_SUCCESS;
}

int Tie::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_resistance = trial_stiffness * t_disp;
    trial_resistance(t_span) -= pseudo_load * penalty * weight_pool;

    return SUANPAN_SUCCESS;
}

int Tie::clear_status() {
    trial_resistance.zeros(get_total_number());
    trial_resistance(t_span) -= pseudo_load * penalty * weight_pool;

    current_resistance = trial_resistance;
    return SUANPAN_SUCCESS;
}

int Tie::commit_status() {
    current_resistance = trial_resistance;
    return SUANPAN_SUCCESS;
}

int Tie::reset_status() {
    trial_resistance = current_resistance;
    return SUANPAN_SUCCESS;
}
