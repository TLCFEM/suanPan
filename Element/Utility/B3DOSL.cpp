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

#include "B3DOSL.h"
#include <Element/Element.h>

const span B3DOSL::sa(0, 2), B3DOSL::sb(3, 5), B3DOSL::sc(7, 9), B3DOSL::sd(10, 12);

OrientationType B3DOSL::get_orientation_type() const { return OrientationType::B3DOS; }

unique_ptr<Orientation> B3DOSL::get_copy() { return make_unique<B3DOSL>(*this); }

vec B3DOSL::to_local_vec(const vec& g_disp) const {
    vec t_disp = g_disp;
    t_disp(sa) = direction_cosine * t_disp(sa);
    t_disp(sb) = direction_cosine * t_disp(sb);
    t_disp(sc) = direction_cosine * t_disp(sc);
    t_disp(sd) = direction_cosine * t_disp(sd);

    vec l_disp(9);

    l_disp(0) = t_disp(7) - t_disp(0); // axial

    auto t_factor = (t_disp(1) - t_disp(8)) / length; // strong axis bending
    l_disp(1) = t_disp(5) + t_factor;                 // near node strong axis
    l_disp(2) = t_disp(12) + t_factor;                // far node strong axis

    t_factor = (t_disp(9) - t_disp(2)) / length; // weak axis bending
    l_disp(3) = t_disp(4) + t_factor;            // near node weak axis
    l_disp(4) = t_disp(11) + t_factor;           // far node weak axis

    l_disp(5) = t_disp(3);  // near node torsion
    l_disp(6) = t_disp(10); // far node torsion

    l_disp(7) = t_disp(6);  // near node warping
    l_disp(8) = t_disp(13); // far node warping

    return l_disp;
}

vec B3DOSL::to_global_vec(const vec& l_disp) const {
    vec g_disp(14, fill::none);

    g_disp(0) = -(g_disp(7) = l_disp(0)); // axial

    g_disp(8) = -(g_disp(1) = ((g_disp(5) = l_disp(1)) + (g_disp(12) = l_disp(2))) / length); // strong axis bending
    g_disp(2) = -(g_disp(9) = ((g_disp(4) = l_disp(3)) + (g_disp(11) = l_disp(4))) / length); // weak axis bending

    g_disp(3) = l_disp(5);
    g_disp(10) = l_disp(6);

    g_disp(6) = l_disp(7);
    g_disp(13) = l_disp(8);

    g_disp(sa) = direction_cosine.t() * g_disp(sa);
    g_disp(sb) = direction_cosine.t() * g_disp(sb);
    g_disp(sc) = direction_cosine.t() * g_disp(sc);
    g_disp(sd) = direction_cosine.t() * g_disp(sd);

    return g_disp;
}

mat B3DOSL::to_global_mass_mat(const mat&) const { return {14, 14, fill::zeros}; }

mat B3DOSL::to_global_stiffness_mat(const mat& l_mat) const {
    mat t_mat(9, 14, fill::zeros); // eq. 2.11

    const auto t_factor = 1. / length;

    // axial
    t_mat(0, 0) = -(t_mat(0, 7) = 1.);
    // strong axis bending near node
    t_mat(1, 8) = -(t_mat(1, 1) = t_factor);
    t_mat(1, 5) = 1.;
    // strong axis bending far node
    t_mat(2, 8) = -(t_mat(2, 1) = t_factor);
    t_mat(2, 12) = 1.;
    // weak axis bending near node
    t_mat(3, 2) = -(t_mat(3, 9) = t_factor);
    t_mat(3, 4) = 1.;
    // weak axis bending far node
    t_mat(4, 2) = -(t_mat(4, 9) = t_factor);
    t_mat(4, 11) = 1.;
    // torsion
    t_mat(5, 3) = t_mat(6, 10) = 1.;
    // warping
    t_mat(7, 6) = t_mat(8, 13) = 1.;

    mat g_mat(14, 14, fill::eye);
    g_mat(sa, sa) = direction_cosine;
    g_mat(sb, sb) = direction_cosine;
    g_mat(sc, sc) = direction_cosine;
    g_mat(sd, sd) = direction_cosine;

    t_mat *= g_mat;

    return t_mat.t() * l_mat * t_mat;
}
