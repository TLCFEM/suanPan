/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "Tag.h"
#include <Toolbox/utility.h>

Tag::Tag(const unsigned T)
    : unique_tag(T) {}

void Tag::set_tag(const unsigned T) const { suanpan::hacker(unique_tag) = T; }

unsigned Tag::get_tag() const { return unique_tag; }

void Tag::enable() { alive = true; }

void Tag::disable() { alive = false; }

void Tag::guard() { guarded = true; }

void Tag::unguard() { guarded = false; }

bool Tag::is_active() const { return alive; }

bool Tag::is_guarded() const { return guarded; }

void Tag::print() {
    suanpan_info("A tagged object.\n");
}
