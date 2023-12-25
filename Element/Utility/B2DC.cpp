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

#include "B2DC.h"
#include <Element/Element.h>
#include <Toolbox/tensor.h>

void B2DC::update_transformation() {
    const auto coord = get_coordinate(element_ptr, 2);
    const auto t_disp = get_trial_displacement(element_ptr);

    vec x_axis(2);
    x_axis(0) = coord(1, 0) - coord(0, 0) + t_disp(3) - t_disp(0);
    x_axis(1) = coord(1, 1) - coord(0, 1) + t_disp(4) - t_disp(1);

    length = norm(x_axis);
    direction_cosine = x_axis / length;
    inclination = transform::atan2(direction_cosine);

    if(original_position.is_empty()) original_position = {length, inclination};

    form_trans_mat(direction_cosine);
}

bool B2DC::is_nlgeom() const { return true; }

unique_ptr<Orientation> B2DC::get_copy() { return make_unique<B2DC>(*this); }

vec B2DC::to_local_vec(const vec& g_vec) const {
    const auto& initial_length = original_position(0);
    const auto& initial_inclination = original_position(1);
    const auto& new_cos = direction_cosine(0);
    const auto& new_sin = direction_cosine(1);

    vec l_vec(3);

    l_vec(0) = length - initial_length;

    auto t_angle = initial_inclination + g_vec(2);
    auto t_sin = sin(t_angle), t_cos = cos(t_angle);

    l_vec(1) = atan((new_cos * t_sin - new_sin * t_cos) / (new_cos * t_cos + new_sin * t_sin));

    t_angle = initial_inclination + g_vec(5);
    t_sin = sin(t_angle);
    t_cos = cos(t_angle);

    l_vec(2) = atan((new_cos * t_sin - new_sin * t_cos) / (new_cos * t_cos + new_sin * t_sin));

    return l_vec;
}

mat B2DC::to_global_geometry_mat(const mat& l_mat) const {
    mat g_mat(6, 6, fill::zeros);

    const auto C2 = direction_cosine(0) * direction_cosine(0) / length;
    const auto S2 = direction_cosine(1) * direction_cosine(1) / length;
    const auto CS = direction_cosine(0) * direction_cosine(1) / length;

    const auto TT = (l_mat(1) + l_mat(2)) / length;
    const auto TA = C2 * l_mat(0) + 2. * CS * TT;
    const auto TB = S2 * l_mat(0) - 2. * CS * TT;
    const auto TC = CS * l_mat(0) + (S2 - C2) * TT;

    g_mat(1, 4) = g_mat(4, 1) = -(g_mat(1, 1) = g_mat(4, 4) = TA);
    g_mat(0, 3) = g_mat(3, 0) = -(g_mat(0, 0) = g_mat(3, 3) = TB);
    g_mat(0, 4) = g_mat(1, 3) = g_mat(3, 1) = g_mat(4, 0) = TC;
    g_mat(0, 1) = g_mat(1, 0) = g_mat(3, 4) = g_mat(4, 3) = -TC;

    return g_mat;
}
