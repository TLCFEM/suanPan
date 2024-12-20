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

#include "PeakOriented.h"

double PeakOriented::compute_compression_residual(const double strain, const double stress) const { return strain - stress / initial_stiffness(0); }

double PeakOriented::compute_tension_residual(const double strain, const double stress) const { return strain - stress / initial_stiffness(0); }

PeakOriented::PeakOriented(const int T, const double R)
    : SimpleHysteresis(T, 1., R) {}

void PeakOriented::print() {
    suanpan_info("A peak oriented hysteresis model.\n");
    SimpleHysteresis::print();
}
