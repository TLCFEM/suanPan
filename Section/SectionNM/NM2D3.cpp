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

#include "NM2D3.h"

double NM2D3::compute_f(const vec& s, const vec& h) const { return compute_sf(s, h); }

vec NM2D3::compute_df(const vec& s, const vec& h) const { return compute_dsf(s, h); }

mat NM2D3::compute_ddf(const vec& s, const vec& h) const { return compute_ddsf(s, h); }

NM2D3::NM2D3(const unsigned T, const double EEA, const double EEIS, const double NP, const double MSP, const double CC, const double HH, const double HS, const double HD, const double KK, const double KB, const double LD, mat&& PS)
    : SurfaceNM2D(CC, std::forward<mat>(PS))
    , VAFNM(T, EEA, EEIS, HH, HS, HD, KK, KB, LD, vec{NP, MSP}) {}

unique_ptr<Section> NM2D3::get_copy() { return make_unique<NM2D3>(*this); }

void NM2D3::print() {
    suanpan_info("A N-M based section. doi: 10.1007/978-94-007-6573-3_3\n");
    current_deformation.t().print("Local Deformation:");
    current_resistance.t().print("Local Resistance:");
}