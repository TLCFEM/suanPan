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

#include "T3DL.h"
#include <Element/Element.h>

const span T3DL::IS(0, 2);
const span T3DL::JS(3, 5);

T3DL::T3DL(const unsigned T)
    : Orientation(T) {}

unique_ptr<Orientation> T3DL::get_copy() { return make_unique<T3DL>(*this); }

void T3DL::update_transformation() {
    if(!direction_cosine.is_empty()) return;

    const auto coord = get_coordinate(element_ptr, 3);

    vec x_axis(3);
    x_axis(0) = coord(1, 0) - coord(0, 0);
    x_axis(1) = coord(1, 1) - coord(0, 1);
    x_axis(2) = coord(1, 2) - coord(0, 2);

    length = norm(x_axis);

    direction_cosine = x_axis / length;
}

vec T3DL::to_local_vec(const vec& g_disp) const { return vec{dot(direction_cosine, g_disp(JS) - g_disp(IS))}; }

vec T3DL::to_global_vec(const vec& l_disp) const {
    vec g_vec(6);

    g_vec(JS) = l_disp(0) * direction_cosine;
    g_vec(IS) = -g_vec(JS);

    return g_vec;
}

mat T3DL::to_global_mass_mat(const mat& l_mat) const {
    mat g_mat(6, 6, fill::zeros);
    g_mat.diag().fill(.5 * length * l_mat(0));
    return g_mat;
}

mat T3DL::to_global_stiffness_mat(const mat& l_mat) const {
    mat g_mat(6, 6);

    g_mat(IS, IS) = l_mat(0) * direction_cosine * direction_cosine.t();
    g_mat(JS, JS) = g_mat(IS, IS);
    g_mat(IS, JS) = -g_mat(IS, IS);
    g_mat(JS, IS) = g_mat(IS, JS);

    return g_mat;
}
