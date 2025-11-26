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

#include "ConditionalModifier.h"

#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Load/Amplitude/Ramp.h>
#include <Step/Step.h>

uvec ConditionalModifier::update_active_dof(const shared_ptr<DomainBase>& D) {
    std::vector<uword> active_dof;

    const auto check = [&](const shared_ptr<Node>& node) {
        if(!node || !node->is_active()) return;
        auto& t_dof = node->get_reordered_dof();
        for(const auto J : dof_reference)
            if(J < t_dof.n_elem) active_dof.emplace_back(t_dof(J));
    };

    if(target_node.is_empty())
        for(auto& node : D->get_node_pool()) check(node);
    else
        for(const auto I : target_node) check(D->get<Node>(I));

    return active_dof;

    // const auto check = [&](const shared_ptr<Node>& node) {
    //     if(!node || !node->is_active()) return;
    //     auto& t_dof = node->get_reordered_dof();
    //     auto& t_identifier = node->get_dof_identifier();
    //     for(auto J = 0u; J < node->get_dof_number(); ++J)
    //         if(dof_identifier.contains(t_identifier[J])) active_dof.emplace_back(t_dof[J]);
    // };
    //
    // if(target_object.is_empty())
    //     for(auto& node : D->get_node_pool()) check(node);
    // else
    //     for(const auto tag : target_object) check(D->get<Node>(tag));
}

double ConditionalModifier::get_amplitude(const shared_ptr<DomainBase>& D) const { return amplitude->get_amplitude(D->get_factory()->get_trial_time()); }

ConditionalModifier::ConditionalModifier(const unsigned T, const unsigned AT, uvec&& N, uvec&& D)
    : UniqueTag(T)
    , amplitude_tag(AT)
    , dof_reference(D - 1)
    , target_node(std::move(N)) {}

int ConditionalModifier::initialize(const shared_ptr<DomainBase>& D) {
    amplitude = D->get<Amplitude>(amplitude_tag);
    if(nullptr == amplitude || !amplitude->is_active()) amplitude = Ramp(0);

    auto start_time = 0.;
    // ReSharper disable once CppUseElementsView
    for(auto& [t_tag, t_step] : D->get_step_pool()) {
        if(t_step->get_tag() >= start_step) break;
        start_time += t_step->get_time_period();
    }
    amplitude->set_start_time(start_time);

    target_dof = update_active_dof(D);

    initialized = true;

    return SUANPAN_SUCCESS;
}

int ConditionalModifier::process_resistance(const shared_ptr<DomainBase>& D) { return process(D); }

const uvec& ConditionalModifier::get_node_encoding() const { return target_node; }

const uvec& ConditionalModifier::get_dof_encoding() const { return target_dof; }

void ConditionalModifier::deinitialize() { initialized = false; }

bool ConditionalModifier::is_initialized() const { return initialized; }

void ConditionalModifier::set_start_step(const unsigned ST) {
    start_step = std::max(1u, ST);
    if(end_step <= start_step) end_step = start_step + 1;
}

unsigned ConditionalModifier::get_start_step() const { return start_step; }

void ConditionalModifier::set_end_step(const unsigned ST) { end_step = ST; }

unsigned ConditionalModifier::get_end_step() const { return end_step; }

bool ConditionalModifier::validate_step(const shared_ptr<DomainBase>& D) const {
    const auto t_step = D->get_current_step_tag();
    return t_step >= start_step && t_step < end_step && is_active();
}

GroupModifier::GroupModifier(uvec&& N)
    : groups(std::move(N)) {}

uvec GroupModifier::update_object_tag(const shared_ptr<DomainBase>& D) const { return D->flatten_group(groups); }
