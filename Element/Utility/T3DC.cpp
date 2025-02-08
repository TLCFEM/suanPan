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

#include "T3DC.h"
#include <Element/Element.h>

void T3DC::update_transformation() {
    const auto coord = get_coordinate(element_ptr, 3);
    const auto t_disp = get_trial_displacement(element_ptr);

    vec x_axis(3);
    x_axis(0) = coord(1, 0) - coord(0, 0) + t_disp(3) - t_disp(0);
    x_axis(1) = coord(1, 1) - coord(0, 1) + t_disp(4) - t_disp(1);
    x_axis(2) = coord(1, 2) - coord(0, 2) + t_disp(5) - t_disp(2);

    length = norm(x_axis);

    direction_cosine = x_axis / length;
}

bool T3DC::is_nlgeom() const { return true; }

unique_ptr<Orientation> T3DC::get_copy() { return make_unique<T3DC>(*this); }

mat T3DC::to_global_geometry_mat(const mat& l_mat) const {
    auto g_mat = to_global_stiffness_mat(-l_mat);

    g_mat.diag() += l_mat(0);
    for(unsigned I = 0; I < 3; ++I) g_mat(I + 3llu, I) = g_mat(I, I + 3llu) -= l_mat(0);

    return g_mat;
}
