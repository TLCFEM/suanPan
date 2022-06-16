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

#include "BodyForce.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Element/Element.h>
#include <Load/Amplitude/Amplitude.h>

BodyForce::BodyForce(const unsigned T, const unsigned S, const double L, uvec&& N, const unsigned D, const unsigned AT)
    : Load(T, S, AT, std::forward<uvec>(N), uvec{D}, L) {}

BodyForce::BodyForce(const unsigned T, const unsigned S, const double L, uvec&& N, uvec&& D, const unsigned AT)
    : Load(T, S, AT, std::forward<uvec>(N), std::forward<uvec>(D), L) {}

int BodyForce::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    trial_load.zeros(W->get_size());

    const auto final_load = pattern * magnitude->get_amplitude(W->get_trial_time());

    for(const auto I : node_encoding)
        if(auto& t_element = D->get<Element>(I); nullptr != t_element && t_element->is_active()) {
            vec t_body_load(t_element->get_dof_number(), fill::zeros);
            for(const auto J : dof_reference) if(J < t_element->get_dof_number()) t_body_load(J) = final_load;
            if(const auto& t_body_force = t_element->update_body_force(t_body_load); !t_body_force.empty()) trial_load(t_element->get_dof_encoding()) += t_body_force;
        }

    return SUANPAN_SUCCESS;
}
