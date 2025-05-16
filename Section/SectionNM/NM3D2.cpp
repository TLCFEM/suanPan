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

#include "NM3D2.h"

double NM3D2::compute_f(const vec& s, const vec& h) const { return compute_sf(s, h); }

vec NM3D2::compute_df(const vec& s, const vec& h) const { return compute_dsf(s, h); }

mat NM3D2::compute_ddf(const vec& s, const vec& h) const { return compute_ddsf(s, h); }

NM3D2::NM3D2(const unsigned T, const double EEA, const double EEIS, const double EEIW, const double NP, const double MSP, const double MWP, const double CC, const double HH, const double KK, const double LD, mat&& PS)
    : SurfaceNM3D(CC, std::move(PS))
    , LinearHardeningNM(T, EEA, EEIS, EEIW, HH, KK, LD, vec{NP, MSP, MWP}) {}

unique_ptr<Section> NM3D2::get_copy() { return std::make_unique<NM3D2>(*this); }
