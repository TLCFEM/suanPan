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

#include "Constraint.h"
#include <Domain/DomainBase.h>
#include <Load/Amplitude/Ramp.h>
#include <Step/Step.h>

const double Constraint::multiplier = 1E8;

Constraint::Constraint(const unsigned T, const unsigned ST, const unsigned AT, uvec&& N, uvec&& D, const unsigned S)
	: Tag(T)
	, start_step(0 == ST ? 1 : ST)
	, amplitude_tag(AT)
	, num_size(S)
	, nodes(std::forward<uvec>(N))
	, dofs(std::forward<uvec>(D)) { suanpan_debug("Constraint %u ctor() called.\n", get_tag()); }

/**
 * \brief default destructor.
 */
Constraint::~Constraint() { suanpan_debug("Constraint %u dtor() called.\n", get_tag()); }

int Constraint::initialize(const shared_ptr<DomainBase>& D) {
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

	set_initialized(true);

	return SUANPAN_SUCCESS;
}

void Constraint::set_initialized(const bool B) const { access::rw(initialized) = B; }

bool Constraint::is_initialized() const { return initialized; }

void Constraint::set_multiplier_size(const unsigned S) { num_size = S; }

unsigned Constraint::get_multiplier_size() const { return num_size; }

/**
 * \brief method to set `start_step`.
 * \param ST `start_step`
 */
void Constraint::set_start_step(const unsigned ST) {
	start_step = ST;
	if(end_step <= start_step) end_step = start_step + 1;
}

/**
 * \brief method to get `start_step`.
 * \return `start_step`
 */
unsigned Constraint::get_start_step() const { return start_step; }

/**
 * \brief method to set `start_step`.
 * \param ST `end_step`
 */
void Constraint::set_end_step(const unsigned ST) { end_step = ST; }

/**
 * \brief method to get `start_step`.
 * \return `end_step`
 */
unsigned Constraint::get_end_step() const { return end_step; }

bool Constraint::validate_step(const shared_ptr<DomainBase>& D) const {
	if(const auto& t_step = D->get_current_step_tag(); t_step < start_step || t_step >= end_step) return false;

	return true;
}

void Constraint::update_incre_lambda(const vec& i_lambda) { trial_lambda += i_lambda; }

void Constraint::commit_status() { current_lambda = trial_lambda; }

void Constraint::clear_status() {
	current_lambda = trial_lambda = 0.;
	access::rw(initialized) = false;
}

void Constraint::reset_status() { trial_lambda = current_lambda; }

void set_constraint_multiplier(const double M) { access::rw(Constraint::multiplier) = M; }
