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

#include "Decay.h"

Decay::Decay(const unsigned T, const double AA, const double TTD)
    : Amplitude(T)
    , A(AA)
    , TD(TTD) {}

unique_ptr<Amplitude> Decay::get_copy() { return std::make_unique<Decay>(*this); }

double Decay::get_amplitude(const double T) {
    const auto step_time = T - start_time;

    return step_time <= 0. ? 0. : A * exp(-start_time / TD);
}

void Decay::print() {
    suanpan_info("Decay.\n");
}
