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

#include "PolyJ2.h"

double PolyJ2::compute_k(const double p_strain) const {
    vec t_vec(poly_para.n_elem);

    t_vec(0) = 1.;
    for(uword I = 1; I < t_vec.n_elem; ++I) t_vec(I) = t_vec(I - 1) * p_strain;

    return dot(poly_para, t_vec);
}

double PolyJ2::compute_dk(const double p_strain) const {
    vec t_vec(poly_para.n_elem);

    t_vec(0) = 0.;
    t_vec(1) = 1.;

    auto t_strain = 1.;
    for(uword I = 2; I < t_vec.n_elem; ++I) t_vec(I) = static_cast<double>(I) * (t_strain *= p_strain);

    return dot(poly_para, t_vec);
}

double PolyJ2::compute_h(const double) const { return 0.; }

double PolyJ2::compute_dh(const double) const { return 0.; }

PolyJ2::PolyJ2(const unsigned T, const double E, const double V, vec&& H, const double R)
    : DataPolyJ2{std::forward<vec>(H)}
    , NonlinearJ2(T, E, V, R) {}

unique_ptr<Material> PolyJ2::get_copy() { return make_unique<PolyJ2>(*this); }

void PolyJ2::print() {
    suanpan_info("A 3D polynomial hardening model.\n");
}
