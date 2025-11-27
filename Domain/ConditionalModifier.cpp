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

#include <Domain/Node.h>
#include <Load/Amplitude/Ramp.h>
#include <Step/Step.h>

bool ConditionalModifier::validate_dof(const shared_ptr<DomainBase>& D) {
    const auto check = [&](const shared_ptr<Node>& node) { return !node || !node->is_active() || !node->validate_dof(dof_order); };

    if(target_node.is_empty())
        for(auto& node : D->get_node_pool()) {
            if(check(node)) return false;
        }
    else
        for(const auto tag : target_node) {
            if(check(D->get<Node>(tag))) return false;
        }

    return true;
}

uvec ConditionalModifier::update_active_dof(const shared_ptr<DomainBase>& D) {
    const auto& ref_component = get_dof_component();

    if(ref_component.empty()) return {};

    std::vector<uword> active_dof;

    const auto check = [&](const shared_ptr<Node>& node) {
        if(!node || !node->is_active()) return;
        suanpan::append_to(active_dof, node->get_dof(ref_component));
    };

    if(target_node.is_empty())
        for(auto& node : D->get_node_pool()) check(node);
    else
        for(const auto tag : target_node) check(D->get<Node>(tag));

    return active_dof;
}

double ConditionalModifier::get_amplitude(const shared_ptr<DomainBase>& D) const { return amplitude->get_amplitude(D->get_factory()->get_trial_time()); }

const std::vector<Node::DOF>& ConditionalModifier::get_dof_component() const { return dof_component.empty() ? dof_order : dof_component; }

ConditionalModifier::ConditionalModifier(const unsigned T, const unsigned AT, uvec&& OT, std::vector<Node::DOF>&& DO, std::vector<Node::DOF>&& DC)
    : UniqueTag(T)
    , amplitude_tag(AT)
    , dof_component(std::move(DC))
    , dof_order(std::move(DO))
    , target_node(std::move(OT)) {}

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

    if(!validate_dof(D)) return SUANPAN_FAIL;

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

std::vector<Node::DOF> parse_dof(const std::string_view token) {
    if(is_equal_any(token, "PINNED", "P")) return std::vector{Node::DOF::U1, Node::DOF::U2, Node::DOF::U3};
    if(is_equal_any(token, "ENCASTRE", "E")) return std::vector{Node::DOF::U1, Node::DOF::U2, Node::DOF::U3, Node::DOF::UR1, Node::DOF::UR2, Node::DOF::UR3};
    if(is_equal_any(token, "XSYMM", "X")) return std::vector{Node::DOF::U1, Node::DOF::UR2, Node::DOF::UR3};
    if(is_equal_any(token, "YSYMM", "Y")) return std::vector{Node::DOF::UR1, Node::DOF::U2, Node::DOF::UR3};
    if(is_equal_any(token, "ZSYMM", "Z")) return std::vector{Node::DOF::UR1, Node::DOF::UR2, Node::DOF::U3};
    if(is_equal_any(token, "1", "U1")) return std::vector{Node::DOF::U1};
    if(is_equal_any(token, "2", "U2")) return std::vector{Node::DOF::U2};
    if(is_equal_any(token, "3", "U3")) return std::vector{Node::DOF::U3};
    if(is_equal_any(token, "4", "U4", "UR1")) return std::vector{Node::DOF::UR1};
    if(is_equal_any(token, "5", "U5", "UR2")) return std::vector{Node::DOF::UR2};
    if(is_equal_any(token, "6", "U6", "UR3")) return std::vector{Node::DOF::UR3};
    if(is_equal(token, "FU1")) return std::vector{Node::DOF::FU1};
    if(is_equal(token, "FU2")) return std::vector{Node::DOF::FU2};
    if(is_equal(token, "FU3")) return std::vector{Node::DOF::FU3};
    if(is_equal(token, "FUR1")) return std::vector{Node::DOF::FUR1};
    if(is_equal(token, "FUR2")) return std::vector{Node::DOF::FUR2};
    if(is_equal(token, "FUR3")) return std::vector{Node::DOF::FUR3};
    if(is_equal(token, "RADIAL")) return std::vector{Node::DOF::RADIAL};
    if(is_equal(token, "AXIAL")) return std::vector{Node::DOF::AXIAL};
    if(is_equal(token, "RS")) return std::vector{Node::DOF::RS};
    if(is_equal(token, "RW")) return std::vector{Node::DOF::RW};
    if(is_equal(token, "DAMAGE")) return std::vector{Node::DOF::DAMAGE};
    if(is_equal(token, "PRESSURE")) return std::vector{Node::DOF::PRESSURE};
    if(is_equal(token, "TEMPERATURE")) return std::vector{Node::DOF::TEMPERATURE};
    if(is_equal(token, "WARP")) return std::vector{Node::DOF::WARP};

    return {};
}
