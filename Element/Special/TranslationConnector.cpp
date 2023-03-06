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

#include "TranslationConnector.h"
#include <Domain/DOF.h>
#include <Domain/DomainBase.h>

TranslationConnector::TranslationConnector(const unsigned T, uvec&& N, const unsigned D, const double P)
    : Element(T, c_node, D, std::forward<uvec>(N), 2u == D ? std::vector{DOF::U1, DOF::U2} : std::vector{DOF::U1, DOF::U2, DOF::U3})
    , c_dof(D)
    , sa(2u == c_dof ? span(0, 1) : span(0, 2))
    , sb(2u == c_dof ? span(2, 3) : span(3, 5))
    , sc(2u == c_dof ? span(4, 5) : span(6, 8))
    , alpha(fabs(P)) {}

int TranslationConnector::initialize(const shared_ptr<DomainBase>&) {
    const mat coor = get_coordinate(c_dof).t();
    const vec x2x1 = coor.col(1) - coor.col(0);
    const vec x3x1 = coor.col(2) - coor.col(0);

    access::rw(s_near) = dot(x3x1, x2x1) / dot(x2x1, x2x1);
    access::rw(s_far) = 1. - s_near;

    mat dvdx(c_dof, c_size, fill::none);
    dvdx.cols(sa) = (alpha - alpha * s_near) * eye(c_dof, c_dof);
    dvdx.cols(sb) = alpha * s_near * eye(c_dof, c_dof);
    dvdx.cols(sc) = -alpha * eye(c_dof, c_dof);

    initial_stiffness.set_size(c_size, c_size);
    initial_stiffness.rows(sa) = -s_far * dvdx;
    initial_stiffness.rows(sb) = -s_near * dvdx;
    initial_stiffness.rows(sc) = dvdx;

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int TranslationConnector::update_status() {
    const mat coor = get_coordinate(c_dof).t();
    const auto t_disp = get_trial_displacement();

    const vec x1 = coor.col(0) + t_disp(sa);
    const vec x2x1 = coor.col(1) + t_disp(sb) - x1;
    const vec x3x1 = coor.col(2) + t_disp(sc) - x1;

    const vec diff_vec = alpha * (x2x1 * s_near - x3x1);

    trial_resistance.set_size(c_size);
    trial_resistance.rows(sa) = -s_far * diff_vec;
    trial_resistance.rows(sb) = -s_near * diff_vec;
    trial_resistance.rows(sc) = diff_vec;

    return SUANPAN_SUCCESS;
}
