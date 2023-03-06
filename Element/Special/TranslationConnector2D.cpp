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

#include "TranslationConnector2D.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/DOF.h>

TranslationConnector2D::TranslationConnector2D(const unsigned T, uvec&& N, const double P)
    : Element(T, c_node, c_dof, std::forward<uvec>(N), {DOF::U1, DOF::U2})
    , alpha(fabs(P)) {}

int TranslationConnector2D::initialize(const shared_ptr<DomainBase>&) {
    const mat coor = get_coordinate(2).t();
    const vec2 x2x1 = coor.col(1) - coor.col(0);
    const vec2 x3x1 = coor.col(2) - coor.col(0);

    access::rw(s_near) = dot(x3x1, x2x1) / dot(x2x1, x2x1);
    access::rw(s_far) = 1. - s_near;

    mat::fixed<2, 6> dvdx;
    dvdx.cols(0, 1) = (alpha - alpha * s_near) * eye(2, 2);
    dvdx.cols(2, 3) = alpha * s_near * eye(2, 2);
    dvdx.cols(4, 5) = -alpha * eye(2, 2);

    initial_stiffness.set_size(c_size, c_size);
    initial_stiffness.rows(0, 1) = -s_far * dvdx;
    initial_stiffness.rows(2, 3) = -s_near * dvdx;
    initial_stiffness.rows(4, 5) = dvdx;

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int TranslationConnector2D::update_status() {
    const mat coor = get_coordinate(2).t();
    const auto t_disp = get_trial_displacement();

    const vec2 x1 = coor.col(0) + t_disp(span(0, 1));
    const vec2 x2x1 = coor.col(1) + t_disp(span(2, 3)) - x1;
    const vec2 x3x1 = coor.col(2) + t_disp(span(4, 5)) - x1;

    const vec2 diff_vec = alpha * (x2x1 * s_near - x3x1);

    trial_resistance.set_size(c_size);
    trial_resistance.rows(0, 1) = -s_far * diff_vec;
    trial_resistance.rows(2, 3) = -s_near * diff_vec;
    trial_resistance.rows(4, 5) = diff_vec;

    return SUANPAN_SUCCESS;
}

int TranslationConnector2D::clear_status() { return SUANPAN_SUCCESS; }

int TranslationConnector2D::commit_status() { return SUANPAN_SUCCESS; }

int TranslationConnector2D::reset_status() { return SUANPAN_SUCCESS; }
