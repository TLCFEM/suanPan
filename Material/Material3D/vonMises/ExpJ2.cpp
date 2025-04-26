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

#include "ExpJ2.h"

double ExpJ2::compute_k(const double p_strain) const {
    const auto b_strain = -b * p_strain;
    return yield_stress * ((1. + a) * exp(b_strain) - a * exp(2. * b_strain));
}

double ExpJ2::compute_dk(const double p_strain) const {
    const auto b_strain = -b * p_strain;
    return yield_stress * b * (2. * a * exp(2. * b_strain) - (1. + a) * exp(b_strain));
}

double ExpJ2::compute_h(const double) const { return 0.; }

double ExpJ2::compute_dh(const double) const { return 0.; }

ExpJ2::ExpJ2(const unsigned T, const double E, const double V, const double YS, const double PA, const double PB, const double R)
    : DataExpJ2{fabs(YS), PA, PB}
    , NonlinearJ2(T, E, V, R) {}

unique_ptr<Material> ExpJ2::get_copy() { return std::make_unique<ExpJ2>(*this); }

void ExpJ2::print() {
    suanpan_info("A 3D exponential hardening model.\n");
}
