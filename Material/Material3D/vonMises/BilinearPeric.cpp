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

#include "BilinearPeric.h"

double BilinearPeric::compute_k(const double p_strain) const { return hardening_modulus >= 0. || p_strain <= -yield_stress / hardening_modulus ? yield_stress + p_strain * hardening_modulus : 0.; }

double BilinearPeric::compute_dk(const double p_strain) const { return hardening_modulus >= 0. || p_strain <= -yield_stress / hardening_modulus ? hardening_modulus : 0.; }

BilinearPeric::BilinearPeric(const unsigned T, const double E, const double V, const double Y, const double H, const double MU, const double EPS, const double R)
    : DataBilinearPeric{Y, H}
    , NonlinearPeric(T, E, V, MU, EPS, R) {}

unique_ptr<Material> BilinearPeric::get_copy() { return make_unique<BilinearPeric>(*this); }
