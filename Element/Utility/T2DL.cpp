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

#include "T2DL.h"
#include <Element/Element.h>
#include <Toolbox/tensor.h>

const span T2DL::IS(0, 1);
const span T2DL::JS(2, 3);

T2DL::T2DL(const unsigned T)
    : Orientation(T) {}

unsigned T2DL::global_size() const { return 2u; }

unsigned T2DL::local_size() const { return 1u; }

unique_ptr<Orientation> T2DL::get_copy() { return make_unique<T2DL>(*this); }

void T2DL::update_transformation() {
    if(!direction_cosine.is_empty()) return;

    const auto coord = get_coordinate(element_ptr, 2);

    vec x_axis(2);
    x_axis(0) = coord(1, 0) - coord(0, 0);
    x_axis(1) = coord(1, 1) - coord(0, 1);

    length = norm(x_axis);
    direction_cosine = x_axis / length;
    inclination = transform::atan2(direction_cosine);
}

vec T2DL::to_local_vec(const vec& g_disp) const { return vec{dot(direction_cosine, g_disp(JS) - g_disp(IS))}; }

vec T2DL::to_global_vec(const vec& l_disp) const {
    vec g_vec(4);

    g_vec(0) = -(g_vec(2) = l_disp(0) * direction_cosine(0));
    g_vec(1) = -(g_vec(3) = l_disp(0) * direction_cosine(1));

    // g_vec(JS) = l_disp(0) * direction_cosine;
    // g_vec(IS) = -g_vec(JS);

    return g_vec;
}

mat T2DL::to_global_mass_mat(const mat& l_mat) const {
    mat g_mat(4, 4, fill::zeros);
    g_mat.diag().fill(.5 * length * l_mat(0));
    return g_mat;
}

mat T2DL::to_global_stiffness_mat(const mat& l_mat) const {
    mat g_mat(4, 4);

    g_mat(IS, IS) = l_mat(0) * direction_cosine * direction_cosine.t();
    g_mat(JS, JS) = g_mat(IS, IS);
    g_mat(IS, JS) = -g_mat(IS, IS);
    g_mat(JS, IS) = g_mat(IS, JS);

    return g_mat;
}
