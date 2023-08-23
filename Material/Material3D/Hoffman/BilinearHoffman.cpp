/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include "BilinearHoffman.h"

double BilinearHoffman::compute_k(const double p_strain) const { return hardening_modulus >= 0. || p_strain <= -1. / hardening_modulus ? 1. + p_strain * hardening_modulus : 0.; }

double BilinearHoffman::compute_dk(const double p_strain) const { return hardening_modulus >= 0. || p_strain <= -1. / hardening_modulus ? hardening_modulus : 0.; }

BilinearHoffman::BilinearHoffman(const unsigned T, vec&& E, vec&& V, vec&& S, const double H, const double R)
    : DataBilinearHoffman{H}
    , NonlinearHoffman(T, std::forward<vec>(E), std::forward<vec>(V), std::forward<vec>(S), R) {}

unique_ptr<Material> BilinearHoffman::get_copy() { return make_unique<BilinearHoffman>(*this); }
