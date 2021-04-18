/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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

B3DL::B3DL(const unsigned T, const double X, const double Y, const double Z)
	: Orientation(T, vec{X, Y, Z}) {}

B3DL::B3DL(const unsigned T, vec&& XYZ)
	: Orientation(T, std::forward<vec>(XYZ)) {}

unique_ptr<Orientation> B3DL::get_copy() { return make_unique<B3DL>(*this); }

void B3DL::update_transformation() {
	if(!direction_cosine.is_empty()) return;

	const mat coord = get_coordinate(element_ptr, 3).t();

	const vec x_axis = coord.col(1) - coord.col(0);

	length = norm(x_axis);

	direction_cosine.resize(3, 3);
	direction_cosine.col(0) = normalise(x_axis);
	direction_cosine.col(1) = normalise(cross(z_axis, x_axis));
	direction_cosine.col(2) = normalise(z_axis);

	if(std::fabs(dot(direction_cosine.col(0), direction_cosine.col(2))) > 1E-4) suanpan_warning("the local z axis of Element %u is not perpendicular to its cord, please check.\n", element_ptr->get_tag());
}

vec B3DL::to_local_vec(const vec& g_disp) const {
	const vec t_disp = vectorise(direction_cosine * reshape(g_disp, 3, 4));

	vec l_disp(6);

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
	vec g_disp(12);

	g_disp(0) = -(g_disp(6) = l_disp(0));
	g_disp(7) = -(g_disp(1) = ((g_disp(5) = l_disp(1)) + (g_disp(11) = l_disp(2))) / length);
	g_disp(2) = -(g_disp(8) = ((g_disp(4) = l_disp(3)) + (g_disp(10) = l_disp(4))) / length);
	g_disp(3) = -(g_disp(9) = l_disp(5));

	return vectorise(direction_cosine.t() * reshape(g_disp, 3, 4));
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
		const span ISPAN(I3, I3 + 2llu);
		for(auto J = 0; J < 4; ++J) {
			const auto J3 = 3 * J;
			const span JSPAN(J3, J3 + 2llu);
			g_mat(ISPAN, JSPAN) = direction_cosine.t() * g_mat(ISPAN, JSPAN) * direction_cosine;
		}
	}

	return g_mat;
}
