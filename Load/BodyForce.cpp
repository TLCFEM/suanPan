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

#include "BodyForce.h"

#include <Domain/Factory.hpp>
#include <Element/Element.h>

BodyForce::BodyForce(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& D, const unsigned AT)
    : Load(T, AT, {}, {}, std::move(D), L)
    , target_element(std::move(N)) {}

int BodyForce::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    trial_load.zeros(W->get_size());

    const auto final_load = magnitude * get_amplitude(D);

    const auto check = [&](const shared_ptr<Element>& element) {
        if(!element || !element->is_active()) return;

        vec t_body_load(element->get_dof_number(), fill::zeros);
        for(const auto J : get_dof_component())
            if(const auto index = static_cast<std::uint8_t>(J) - 1u; index < element->get_dof_number()) t_body_load(index) = final_load;
        if(auto& t_body_force = element->update_body_force(t_body_load); !t_body_force.empty()) trial_load(element->get_dof_encoding()) += t_body_force; };

    for(const auto tag : target_element) check(D->get<Element>(tag));

    return SUANPAN_SUCCESS;
}

GroupBodyForce::GroupBodyForce(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& D, const unsigned AT)
    : GroupModifier(std::move(N))
    , BodyForce(T, L, {}, std::move(D), AT) {}

int GroupBodyForce::initialize(const shared_ptr<DomainBase>& D) {
    target_element = update_object_tag(D);

    return BodyForce::initialize(D);
}
