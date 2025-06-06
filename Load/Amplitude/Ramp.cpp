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

#include "Ramp.h"

unique_ptr<Amplitude> Ramp::get_copy() { return std::make_unique<Ramp>(*this); }

double Ramp::get_amplitude(const double T) {
    const auto step_time = T - start_time;
    if(step_time < 0.) return 0.;
    if(step_time > 1.) return 1.;
    return step_time;
}

void Ramp::print() {
    suanpan_info("Ramp.\n");
}
