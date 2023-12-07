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

#include "B3DOSC.h"
#include <Element/Element.h>
#include <Toolbox/tensor.h>

B3DOSC::B3DOSC(const unsigned T, vec&& XYZ)
    : B3DC(T, std::move(XYZ)) {
    access::rw(sc) = span(7, 9);
    access::rw(sd) = span(10, 12);
}

OrientationType B3DOSC::get_orientation_type() const { return OrientationType::B3DOS; }

void B3DOSC::update_transformation() {
    const mat t_coor = get_coordinate(element_ptr, 3).t();
    const vec x_axis = t_coor.col(1) - t_coor.col(0);

    update_direct_cosine(x_axis);

    const mat trial_disp = reshape(get_trial_displacement(element_ptr), 7, 2);
    const vec incre_disp = trial_disp.head_rows(3).col(1) - trial_disp.head_rows(3).col(0);
    const vec trial_cord = x_axis + incre_disp;

    // new length
    length = norm(trial_cord);

    elongation = dot(x_axis + trial_cord, incre_disp) / (length + initial_length); // eq. 4.98

    const vec incre_ri = trial_disp.rows(sb).col(0) - trial_rotation.col(0);
    const vec incre_rj = trial_disp.rows(sb).col(1) - trial_rotation.col(1);
    trial_rotation = trial_disp.rows(sb);

    // nodal frame
    trial_n.head_cols(3) = transform::rodrigues(incre_ri) * trial_n.head_cols(3);
    trial_n.tail_cols(3) = transform::rodrigues(incre_rj) * trial_n.tail_cols(3);

    // reference frame
    trial_ref = transform::rodrigues((.5 * (incre_ri + incre_rj)).eval()) * trial_ref;

    // basic deformed frame
    update_e(trial_cord);

    theta.set_size(8);
    theta(6) = trial_disp(6, 0);
    theta(7) = trial_disp(6, 1);
    update_theta();

    const auto a = compute_a();
    const auto lr1 = compute_l(a, r(1)), lr2 = compute_l(a, r(2));
    for(auto I = 0u; I < 6u; ++I) sn(I) = transform::skew_symm(trial_n.col(I));

    const mat na = trial_n.t() * a;

    transformation.zeros(9, 14);

    // row 1
    transformation.row(0).cols(0, 2) = -e(0).t();
    transformation.row(0).cols(7, 9) = e(0).t();

    // row 2
    transformation.row(1) = ni(0).t() * lr1;
    transformation.row(1).cols(0, 2) += na.row(1);
    transformation.row(1).cols(3, 5) += e(0).t() * sni(1) - e(1).t() * sni(0);
    transformation.row(1).cols(7, 9) -= na.row(1);
    transformation.row(1) *= .5 / std::cos(theta(2));

    // row 3
    transformation.row(2) = nj(0).t() * lr1;
    transformation.row(2).cols(0, 2) += na.row(4);
    transformation.row(2).cols(7, 9) -= na.row(4);
    transformation.row(2).cols(10, 12) += e(0).t() * snj(1) - e(1).t() * snj(0);
    transformation.row(2) *= .5 / std::cos(theta(5));

    // row 4
    transformation.row(3) = -ni(0).t() * lr2;
    transformation.row(3).cols(0, 2) -= na.row(2);
    transformation.row(3).cols(3, 5) -= e(0).t() * sni(2) - e(2).t() * sni(0);
    transformation.row(3).cols(7, 9) += na.row(2);
    transformation.row(3) *= .5 / std::cos(theta(1));

    // row 5
    transformation.row(4) = -nj(0).t() * lr2;
    transformation.row(4).cols(0, 2) -= na.row(5);
    transformation.row(4).cols(7, 9) += na.row(5);
    transformation.row(4).cols(10, 12) -= e(0).t() * snj(2) - e(2).t() * snj(0);
    transformation.row(4) *= .5 / std::cos(theta(4));

    // row 6
    transformation.row(5) = ni(1).t() * lr2 - ni(2).t() * lr1;
    transformation.row(5).cols(3, 5) += e(1).t() * sni(2) - e(2).t() * sni(1);
    transformation.row(5) *= .5 / std::cos(theta(0));

    // row 7
    transformation.row(6) = nj(1).t() * lr2 - nj(2).t() * lr1;
    transformation.row(6).cols(10, 12) += e(1).t() * snj(2) - e(2).t() * snj(1);
    transformation.row(6) *= .5 / std::cos(theta(3));

    // rows 8 and 9
    transformation(7, 6) = transformation(8, 13) = 1.;
}

unsigned B3DOSC::nodal_size() const { return 7u; }

unique_ptr<Orientation> B3DOSC::get_copy() { return make_unique<B3DOSC>(*this); }

vec B3DOSC::to_local_vec(const vec&) const { return {elongation, theta(2), theta(5), theta(1), theta(4), theta(0), theta(3), theta(6), theta(7)}; }

mat B3DOSC::to_global_geometry_mat(const mat& l_force) const {
    mat geometry(14, 14, fill::zeros);

    const auto &p1 = l_force(0), &p2 = l_force(1), &p3 = l_force(2), &p4 = l_force(3), &p5 = l_force(4), &p6 = l_force(5), &p7 = l_force(6);

    const auto m2 = p2 * .5 / std::cos(theta(2));
    const auto m3 = p3 * .5 / std::cos(theta(5));
    const auto m4 = p4 * .5 / std::cos(theta(1));
    const auto m5 = p5 * .5 / std::cos(theta(4));
    const auto m6 = p6 * .5 / std::cos(theta(0));
    const auto m7 = p7 * .5 / std::cos(theta(3));

    const mat a = compute_a();

    // KA eq. 4.134
    geometry(sa, sa) = p1 * a;
    geometry(sc, sc) = geometry(sa, sa);
    geometry(sa, sc) = -geometry(sa, sa);
    geometry(sc, sa) = -geometry(sa, sa);

    // KB eq. 4.135
    geometry += p2 * std::tan(theta(2)) * transformation.row(1).t() * transformation.row(1);
    geometry += p3 * std::tan(theta(5)) * transformation.row(2).t() * transformation.row(2);
    geometry += p4 * std::tan(theta(1)) * transformation.row(3).t() * transformation.row(3);
    geometry += p5 * std::tan(theta(4)) * transformation.row(4).t() * transformation.row(4);
    geometry += p6 * std::tan(theta(0)) * transformation.row(5).t() * transformation.row(5);
    geometry += p7 * std::tan(theta(3)) * transformation.row(6).t() * transformation.row(6);

    // KC !!!eq. 4.136 is wrong, m6i and m6j need to be swapped!!!
    geometry += m2 * compute_g(a, r(1), ni(0)) + m3 * compute_g(a, r(1), nj(0)) - m4 * compute_g(a, r(2), ni(0)) - m5 * compute_g(a, r(2), nj(0)) + m6 * (compute_g(a, r(2), ni(1)) - compute_g(a, r(1), ni(2))) + m7 * (compute_g(a, r(2), nj(1)) - compute_g(a, r(1), nj(2)));

    // KD eq. 4.137
    const mat KD2 = compute_l(a, r(2)).t() * (m4 * sni(0) - m6 * sni(1)) - compute_l(a, r(1)).t() * (m2 * sni(0) - m6 * sni(2));
    const mat KD4 = compute_l(a, r(2)).t() * (m5 * snj(0) - m7 * snj(1)) - compute_l(a, r(1)).t() * (m3 * snj(0) - m7 * snj(2));

    geometry.cols(sb) += KD2;
    geometry.cols(sd) += KD4;

    // KE eq. 4.139
    geometry.rows(sb) += KD2.t();
    geometry.rows(sd) += KD4.t();

    // KF eq. 4.140
    const mat KF11 = -m2 * compute_m(a, ni(1)) - m3 * compute_m(a, nj(1)) + m4 * compute_m(a, ni(2)) + m5 * compute_m(a, nj(2));
    geometry(sa, sa) += KF11;
    geometry(sc, sc) += KF11;
    geometry(sa, sc) -= KF11;
    geometry(sc, sa) -= KF11;

    const mat KF12 = -m2 * a * sni(1) + m4 * a * sni(2);
    geometry(sa, sb) += KF12;
    geometry(sb, sa) += KF12.t();
    geometry(sb, sc) -= KF12.t();
    geometry(sc, sb) -= KF12;

    const mat KF14 = -m3 * a * snj(1) + m5 * a * snj(2);
    geometry(sa, sd) += KF14;
    geometry(sd, sa) += KF14.t();
    geometry(sc, sd) -= KF14;
    geometry(sd, sc) -= KF14.t();

    // KF22 eq. 4.141
    geometry(sb, sb) += m2 * (se(1) * sni(0) - se(0) * sni(1)) - m4 * (se(2) * sni(0) - se(0) * sni(2)) + m6 * (se(2) * sni(1) - se(1) * sni(2));
    // KF44 eq. 4.141
    geometry(sd, sd) += m3 * (se(1) * snj(0) - se(0) * snj(1)) - m5 * (se(2) * snj(0) - se(0) * snj(2)) + m7 * (se(2) * snj(1) - se(1) * snj(2));

    return geometry;
}
