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

#include "B2DL.h"
#include <Element/Element.h>
#include <Toolbox/tensor.h>

void B2DL::form_trans_mat(const vec& d_cosine) {
    trans_mat.zeros(3, 6);
    trans_mat(0, 0) = -(trans_mat(0, 3) = d_cosine(0));
    trans_mat(0, 1) = -(trans_mat(0, 4) = d_cosine(1));
    trans_mat(1, 0) = trans_mat(2, 0) = -(trans_mat(1, 3) = trans_mat(2, 3) = d_cosine(1) / length);
    trans_mat(1, 4) = trans_mat(2, 4) = -(trans_mat(1, 1) = trans_mat(2, 1) = d_cosine(0) / length);
    trans_mat(1, 2) = trans_mat(2, 5) = 1.;
}

B2DL::B2DL(const unsigned T, const double X, const double Y, const double Z)
    : Orientation(T, vec{X, Y, Z}) {}

B2DL::B2DL(const unsigned T, vec&& XYZ)
    : Orientation(T, std::forward<vec>(XYZ)) {}

unique_ptr<Orientation> B2DL::get_copy() { return make_unique<B2DL>(*this); }

void B2DL::update_transformation() {
    if(!direction_cosine.is_empty()) return;

    const mat coord = get_coordinate(element_ptr, 2).t();

    const vec x_axis = coord.col(1) - coord.col(0);

    length = norm(x_axis);
    direction_cosine = x_axis / length;
    inclination = transform::atan2(direction_cosine);

    form_trans_mat(direction_cosine);
}

vec B2DL::to_local_vec(const vec& g_disp) const { return trans_mat * g_disp; }

vec B2DL::to_global_vec(const vec& l_disp) const { return trans_mat.t() * l_disp; }

mat B2DL::to_global_mass_mat(const mat& l_mat) const {
    mat g_mat(6, 6, fill::zeros);

    if(l_mat(0) != 0.) {
        mat trans(6, 6, fill::zeros);
        trans(5, 5) = trans(2, 2) = 1.;
        trans(0, 0) = trans(1, 1) = trans(3, 3) = trans(4, 4) = direction_cosine(0);
        trans(1, 0) = trans(4, 3) = -(trans(0, 1) = trans(3, 4) = direction_cosine(1));

        g_mat(1, 1) = g_mat(4, 4) = 156.;
        g_mat(1, 4) = g_mat(4, 1) = 54.;
        g_mat(4, 5) = g_mat(5, 4) = -(g_mat(2, 1) = g_mat(1, 2) = 22. * length);
        g_mat(2, 4) = g_mat(4, 2) = -(g_mat(5, 1) = g_mat(1, 5) = -13. * length);
        g_mat(5, 2) = g_mat(2, 5) = -.75 * (g_mat(5, 5) = g_mat(2, 2) = 4. * length * length);
        g_mat(3, 3) = g_mat(0, 0) = 2. * (g_mat(3, 0) = g_mat(0, 3) = 140.);

        g_mat = trans.t() * g_mat * trans;
    }

    return g_mat *= l_mat(0) * length / 420.;
}

mat B2DL::to_global_stiffness_mat(const mat& l_mat) const { return trans_mat.t() * l_mat * trans_mat; }
