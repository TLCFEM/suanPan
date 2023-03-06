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

void TranslationConnector2D::update_position() {
    const mat coor = get_coordinate(2).t();
    const auto t_disp = get_trial_displacement();
    const vec2 x1 = coor.col(0) + t_disp(span(0, 1));
    const vec2 x2 = coor.col(1) + t_disp(span(2, 3));
    const vec2 x3 = coor.col(2) + t_disp(span(4, 5));
    const vec2 x2x1 = x2 - x1;
    const vec2 x3x1 = x3 - x1;
    const vec2 normalise_x2x1 = normalise(x2x1);
    const auto norm_x2x1 = norm(x2x1);

    const vec2 diff_vec = alpha * (x2x1 * s_near - norm_x2x1 * x3x1);

    mat::fixed<2, 6> dvec;
    dvec.cols(0, 1) = alpha * ((norm_x2x1 - s_near) * eye(2, 2) + x3x1 * normalise_x2x1.t());
    dvec.cols(2, 3) = alpha * (s_near * eye(2, 2) - x3x1 * normalise_x2x1.t());
    dvec.cols(4, 5) = -alpha * norm_x2x1 * eye(2, 2);

    trial_resistance.set_size(c_size);
    trial_resistance.rows(0, 1) = -s_far * diff_vec;
    trial_resistance.rows(2, 3) = -s_near * diff_vec;
    trial_resistance.rows(4, 5) = diff_vec;

    trial_stiffness.set_size(c_size, c_size);
    trial_stiffness.rows(0, 1) = -s_far * dvec;
    trial_stiffness.rows(2, 3) = -s_near * dvec;
    trial_stiffness.rows(4, 5) = dvec;
}

TranslationConnector2D::TranslationConnector2D(const unsigned T, uvec&& N, const double P)
    : Element(T, c_node, c_dof, std::forward<uvec>(N), {DOF::U1, DOF::U2})
    , alpha(P) {}

int TranslationConnector2D::initialize(const shared_ptr<DomainBase>&) {
    initial_stiffness.zeros(c_size, c_size);

    const mat coor = get_coordinate(2).t();
    const vec2 x1 = coor.col(0);
    const vec2 x2x1 = coor.col(1) - x1;
    const vec2 x3x1 = coor.col(2) - x1;

    access::rw(s_near) = dot(x3x1, x2x1) / dot(x2x1, x2x1);
    access::rw(s_far) = 1. - s_near;

    update_position();

    current_stiffness = initial_stiffness = trial_stiffness;

    return SUANPAN_SUCCESS;
}

int TranslationConnector2D::update_status() {
    update_position();

    return SUANPAN_SUCCESS;
}

int TranslationConnector2D::clear_status() { return SUANPAN_SUCCESS; }

int TranslationConnector2D::commit_status() { return SUANPAN_SUCCESS; }

int TranslationConnector2D::reset_status() { return SUANPAN_SUCCESS; }
