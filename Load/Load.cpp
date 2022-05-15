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

#include "Load.h"

constexpr double Load::multiplier = 1E8;

Load::Load(const unsigned T, const unsigned ST, const unsigned AT, uvec&& NT, uvec&& DT, const double PT)
    : ConditionalModifier(T, ST, AT, std::forward<uvec>(NT), std::forward<uvec>(DT))
    , pattern(PT) { suanpan_debug("Load %u ctor() called.\n", get_tag()); }

Load::~Load() { suanpan_debug("Load %u dtor() called.\n", get_tag()); }

void Load::enable_displacement_control() const { access::rw(mpdc_flag) = true; }

bool Load::if_displacement_control() const { return mpdc_flag; }

const vec& Load::get_trial_load() const { return trial_load; }

const vec& Load::get_trial_settlement() const { return trial_settlement; }

void set_load_multiplier(const double M) { access::rw(Load::multiplier) = M; }
