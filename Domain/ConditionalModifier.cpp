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

uvec ConditionalModifier::get_nodal_active_dof(const shared_ptr<DomainBase>& D) {
    std::vector<uword> active_dof;
    active_dof.reserve(node_encoding.n_elem * dof_reference.n_elem);

    for(const auto I : node_encoding)
        if(auto& t_node = D->get<Node>(I); nullptr != t_node && t_node->is_active()) {
            auto& t_dof = t_node->get_reordered_dof();
            for(const auto J : dof_reference) if(J < t_dof.n_elem) active_dof.emplace_back(t_dof(J));
        }

    return active_dof;
}

uvec ConditionalModifier::get_all_nodal_active_dof(const shared_ptr<DomainBase>& D) {
    std::vector<uword> active_dof;
    active_dof.reserve(D->get_node() * dof_reference.n_elem);

    for(const auto& I : D->get_node_pool()) {
        auto& t_dof = I->get_reordered_dof();
        for(const auto J : dof_reference) if(J < t_dof.n_elem) active_dof.emplace_back(t_dof(J));
    }

    return active_dof;
}

ConditionalModifier::ConditionalModifier(const unsigned T, const unsigned ST, const unsigned AT, uvec&& N, uvec&& D)
    : UniqueTag(T)
    , start_step(std::max(1u, ST))
    , amplitude_tag(AT)
    , node_encoding(std::move(N))
    , dof_reference(D - 1) {}

int ConditionalModifier::initialize(const shared_ptr<DomainBase>& D) {
    if(0 == amplitude_tag) magnitude = make_shared<Ramp>(0);
    else {
        magnitude = D->get<Amplitude>(amplitude_tag);
        if(nullptr == magnitude || !magnitude->is_active()) magnitude = make_shared<Ramp>(0);
    }

    auto start_time = 0.;
    for(const auto& [t_tag, t_step] : D->get_step_pool()) {
        if(t_step->get_tag() >= start_step) break;
        start_time += t_step->get_time_period();
    }

    magnitude->set_start_step(start_step);
    magnitude->set_start_time(start_time);

    set_initialized(true);

    return SUANPAN_SUCCESS;
}

int ConditionalModifier::process_resistance(const shared_ptr<DomainBase>& D) { return process(D); }

const uvec& ConditionalModifier::get_node_encoding() const { return node_encoding; }

const uvec& ConditionalModifier::get_dof_encoding() const { return dof_encoding; }

void ConditionalModifier::set_initialized(const bool B) const { access::rw(initialized) = B; }

bool ConditionalModifier::is_initialized() const { return initialized; }

void ConditionalModifier::set_start_step(const unsigned ST) {
    start_step = std::max(1u, ST);
    if(end_step <= start_step) end_step = start_step + 1;
}

unsigned ConditionalModifier::get_start_step() const { return start_step; }

void ConditionalModifier::set_end_step(const unsigned ST) { end_step = ST; }

unsigned ConditionalModifier::get_end_step() const { return end_step; }

void ConditionalModifier::set_connected(const bool B) const { access::rw(connected) = B; }

bool ConditionalModifier::is_connected() const { return connected; }

bool ConditionalModifier::validate_step(const shared_ptr<DomainBase>& D) const {
    const auto t_step = D->get_current_step_tag();
    return t_step >= start_step && t_step < end_step && is_active();
}

void ConditionalModifier::update_status(const vec&) {}

void ConditionalModifier::commit_status() {}

void ConditionalModifier::clear_status() { set_initialized(false); }

void ConditionalModifier::reset_status() {}
