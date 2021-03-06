////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2021 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "Load.h"
#include <Domain/DomainBase.h>
#include <Load/Amplitude/Ramp.h>
#include <Step/Step.h>

const double Load::multiplier = 1E8;

Load::Load(const unsigned T, const unsigned ST, const unsigned AT, uvec&& NT, uvec&& DT, const double PT)
	: Tag(T)
	, start_step(0 == ST ? 1 : ST)
	, amplitude_tag(AT)
	, nodes(std::forward<uvec>(NT))
	, dofs(std::forward<uvec>(DT))
	, pattern(PT) { suanpan_debug("Load %u ctor() called.\n", get_tag()); }

Load::~Load() { suanpan_debug("Load %u dtor() called.\n", get_tag()); }

int Load::initialize(const shared_ptr<DomainBase>& D) {
	0 == amplitude_tag ? magnitude = make_shared<Ramp>(0) : magnitude = D->get<Amplitude>(amplitude_tag);

	if(nullptr != magnitude && !magnitude->is_active()) magnitude = nullptr;

	if(nullptr == magnitude) magnitude = make_shared<Amplitude>();

	auto start_time = 0.;
	for(const auto& [t_tag, t_step] : D->get_step_pool()) {
		if(t_step->get_tag() >= start_step) break;
		start_time += t_step->get_time_period();
	}

	magnitude->set_start_step(start_step);
	magnitude->set_start_time(start_time);

	access::rw(initialized) = true;

	return SUANPAN_SUCCESS;
}

void Load::set_initialized(const bool B) const { access::rw(initialized) = B; }

void Load::set_start_step(const unsigned T) {
	start_step = T;
	if(end_step <= start_step) end_step = start_step + 1;
}

unsigned Load::get_start_step() const { return start_step; }

void Load::set_end_step(const unsigned T) { end_step = T; }

unsigned Load::get_end_step() const { return end_step; }

bool Load::validate_step(const shared_ptr<DomainBase>& D) const {
	const auto t_step = D->get_current_step_tag();
	return t_step >= start_step && t_step < end_step;
}

void Load::commit_status() {}

void Load::clear_status() { access::rw(initialized) = false; }

void Load::reset_status() {}

void set_load_multiplier(const double M) { access::rw(Load::multiplier) = M; }
