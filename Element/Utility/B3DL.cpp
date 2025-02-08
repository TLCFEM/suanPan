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

#include "B3DL.h"
#include <Element/Element.h>

OrientationType B3DL::get_orientation_type() const { return OrientationType::B3D; }

unique_ptr<Orientation> B3DL::get_copy() { return make_unique<B3DL>(*this); }

void B3DL::update_transformation() {
    if(!direction_cosine.is_empty()) return;

    const mat coor = get_coordinate(element_ptr, 3).t();

    const vec x_axis = coor.col(1) - coor.col(0);

    length = norm(x_axis);

    direction_cosine.resize(3, 3);
    direction_cosine.col(0) = normalise(x_axis);
    direction_cosine.col(1) = normalise(cross(z_axis, x_axis));
    direction_cosine.col(2) = normalise(cross(x_axis, direction_cosine.col(1)));
}

vec B3DL::to_local_vec(const vec& g_disp) const {
    vec t_disp(g_disp.n_elem, fill::none);
    for(auto I = 0, J = 2; I < 12; I += 3, J += 3) {
        const span sa(I, J);
        t_disp(sa) = direction_cosine * g_disp(sa);
    }

    vec l_disp(6); // eq. 2.11

    l_disp(0) = t_disp(6) - t_disp(0);
    l_disp(5) = t_disp(9) - t_disp(3);
    auto t_factor = (t_disp(1) - t_disp(7)) / length;
    l_disp(1) = t_disp(5) + t_factor;
    l_disp(2) = t_disp(11) + t_factor;
    t_factor = (t_disp(8) - t_disp(2)) / length;
    l_disp(3) = t_disp(4) + t_factor;
    l_disp(4) = t_disp(10) + t_factor;

    return l_disp;
}

vec B3DL::to_global_vec(const vec& l_disp) const {
    vec g_disp(12, fill::none);

    g_disp(0) = -(g_disp(6) = l_disp(0));
    g_disp(7) = -(g_disp(1) = ((g_disp(5) = l_disp(1)) + (g_disp(11) = l_disp(2))) / length);
    g_disp(2) = -(g_disp(8) = ((g_disp(4) = l_disp(3)) + (g_disp(10) = l_disp(4))) / length);
    g_disp(3) = -(g_disp(9) = l_disp(5));

    for(auto I = 0, J = 2; I < 12; I += 3, J += 3) {
        const span sa(I, J);
        g_disp(sa) = direction_cosine.t() * g_disp(sa);
    }

    return g_disp;
}

mat B3DL::to_global_mass_mat(const mat& l_mat) const {
    mat g_mat(12, 12, fill::zeros);

    if(l_mat(0) > 0.) {
        const double fa = l_mat(0) * length / 420.;
        const double fb = fa * length;
        const double fc = fb * length;

        g_mat(0, 0) = g_mat(6, 6) = 140. * fa;
        g_mat(0, 6) = g_mat(6, 0) = 70. * fa;

        g_mat(2, 2) = g_mat(8, 8) = 156. * fa;
        g_mat(2, 8) = g_mat(8, 2) = 54. * fa;
        g_mat(4, 4) = g_mat(10, 10) = 4. * fc;
        g_mat(4, 10) = g_mat(10, 4) = -3. * fc;
        g_mat(2, 4) = g_mat(4, 2) = -22. * fb;
        g_mat(8, 10) = g_mat(10, 8) = -g_mat(2, 4);
        g_mat(2, 10) = g_mat(10, 2) = 13. * fb;
        g_mat(4, 8) = g_mat(8, 4) = -g_mat(2, 10);

        g_mat(1, 1) = g_mat(7, 7) = 156. * fa;
        g_mat(1, 7) = g_mat(7, 1) = 54. * fa;
        g_mat(5, 5) = g_mat(11, 11) = 4. * fc;
        g_mat(5, 11) = g_mat(11, 5) = -3. * fc;
        g_mat(1, 5) = g_mat(5, 1) = 22. * fb;
        g_mat(7, 11) = g_mat(11, 7) = -g_mat(1, 5);
        g_mat(1, 11) = g_mat(11, 1) = -fb * 13.;
        g_mat(5, 7) = g_mat(7, 5) = -g_mat(1, 11);

        for(auto I = 0; I < 4; ++I) {
            const auto I3 = 3 * I;
            const span sa(I3, I3 + 2llu);
            for(auto J = 0; J < 4; ++J) {
                const auto J3 = 3 * J;
                const span sb(J3, J3 + 2llu);
                g_mat(sa, sb) = direction_cosine.t() * g_mat(sa, sb) * direction_cosine;
            }
        }
    }

    return g_mat;
}

mat B3DL::to_global_stiffness_mat(const mat& l_mat) const {
    mat t_mat(6, 12), g_mat(12, 12);

    for(auto I = 0; I < 6; ++I) {
        t_mat(I, 0) = -(t_mat(I, 6) = l_mat(I, 0));
        t_mat(I, 7) = -(t_mat(I, 1) = ((t_mat(I, 5) = l_mat(I, 1)) + (t_mat(I, 11) = l_mat(I, 2))) / length);
        t_mat(I, 2) = -(t_mat(I, 8) = ((t_mat(I, 4) = l_mat(I, 3)) + (t_mat(I, 10) = l_mat(I, 4))) / length);
        t_mat(I, 3) = -(t_mat(I, 9) = l_mat(I, 5));
    }

    for(auto I = 0; I < 12; ++I) {
        g_mat(0, I) = -(g_mat(6, I) = t_mat(0, I));
        g_mat(7, I) = -(g_mat(1, I) = ((g_mat(5, I) = t_mat(1, I)) + (g_mat(11, I) = t_mat(2, I))) / length);
        g_mat(2, I) = -(g_mat(8, I) = ((g_mat(4, I) = t_mat(3, I)) + (g_mat(10, I) = t_mat(4, I))) / length);
        g_mat(3, I) = -(g_mat(9, I) = t_mat(5, I));
    }

    for(auto I = 0; I < 4; ++I) {
        const auto I3 = 3 * I;
        const span sa(I3, I3 + 2llu);
        for(auto J = 0; J < 4; ++J) {
            const auto J3 = 3 * J;
            const span sb(J3, J3 + 2llu);
            g_mat(sa, sb) = direction_cosine.t() * g_mat(sa, sb) * direction_cosine;
        }
    }

    return g_mat;
}
