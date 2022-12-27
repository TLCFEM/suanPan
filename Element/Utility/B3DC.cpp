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

#include "B3DC.h"
#include <Element/Element.h>
#include <Toolbox/tensorToolbox.h>

unique_ptr<Orientation> B3DC::get_copy() { return make_unique<B3DC>(*this); }

void B3DC::commit_status() { current_n = trial_n; }

void B3DC::reset_status() { trial_n = current_n; }

void B3DC::clear_status() {
    direction_cosine.clear();

    trial_n = current_n = direction_cosine;
}

mat B3DC::compute_a() const {
    return (eye(3, 3) - e(0) * e(0).t()) / length; // eq. 4.109
}

mat B3DC::compute_l(const mat& a, const subview_col<double>& rk) const {
    const vec e0r0 = .25 * (e(0) + r(0));
    const mat srk = transform::skew_symm(rk);
    const auto rke0 = .25 * dot(rk, e(0));

    mat l(3, 12, fill::none);

    l.cols(0, 2) = (2. * rke0 * eye(3, 3) + 2. * e0r0 * rk.t()) * a;                             // eq. 4.110
    l.cols(3, 5) = rke0 * transform::skew_symm(r(0)) + (e0r0 * e(0).t() - .5 * eye(3, 3)) * srk; // eq. 4.110

    l.cols(6, 8) = -l.cols(0, 2);
    l.cols(9, 11) = l.cols(3, 5);

    return l;
}

mat B3DC::compute_m(const mat& a, const subview_col<double>& z) const {
    const mat aze0 = a * z * e(0).t();
    return -(aze0 + aze0.t() + a * dot(e(0), z)) / length; // eq. 4.131
}

mat B3DC::compute_g(const mat& a, const subview_col<double>& rk, const subview_col<double>& z) const {
    mat g(12, 12, fill::none);

    const auto srk = transform::skew_symm(rk);
    const auto sr0 = transform::skew_symm(r(0));
    const auto sz = transform::skew_symm(z);
    const auto rke0 = dot(rk, e(0));
    const auto e0r0z = dot(e(0) + r(0), z);
    const mat zrk = z * rk.t();
    const mat ze0 = z * e(0).t();

    const auto sa = span(0, 2), sb = span(3, 5), sc = span(6, 8), sd = span(9, 11);

    g(sa, sa) = -.5 * (a * (zrk + zrk.t()) * a + rke0 * compute_m(a, z) + e0r0z * compute_m(a, rk));
    g(sc, sc) = g(sa, sa);
    g(sa, sc) = -g(sa, sa);
    g(sc, sa) = -g(sa, sa);

    g(sb, sb) = .125 * ((srk * ze0.t() - rke0 * sz) * sr0 + (2. * sz + sr0 * ze0 - e0r0z * se(0)) * srk);
    g(sd, sd) = g(sb, sb);
    g(sb, sd) = g(sb, sb);
    g(sd, sb) = g(sb, sb);

    g(sa, sb) = -.25 * a * ((ze0 + eye(3, 3) * e0r0z) * srk + rk * z.t() * sr0);
    g(sa, sd) = g(sa, sb);
    g(sc, sb) = -g(sa, sb);
    g(sc, sd) = -g(sa, sb);
    g(sb, sa) = g(sa, sb).t();
    g(sb, sc) = -g(sa, sb).t();
    g(sd, sa) = g(sa, sb).t();
    g(sd, sc) = -g(sa, sb).t();

    return g;
}

subview_col<double> B3DC::e(const uword I) const { return basic.col(I); }

subview_col<double> B3DC::r(const uword I) const { return reference.col(I); }

subview_col<double> B3DC::ni(const uword I) const { return trial_n.col(I); }

subview_col<double> B3DC::nj(const uword I) const { return trial_n.col(I + 3); }

const mat& B3DC::sni(const uword I) const { return sn(I); }

const mat& B3DC::snj(const uword I) const { return sn(I + 3); }

void B3DC::update_transformation() {
    const mat t_coor = get_coordinate(element_ptr, 3).t();
    const vec x_axis = t_coor.col(1) - t_coor.col(0);

    if(direction_cosine.empty()) {
        // initial undeformed frame
        direction_cosine.resize(3, 3);
        direction_cosine.col(0) = normalise(x_axis);
        direction_cosine.col(1) = normalise(cross(z_axis, x_axis));
        direction_cosine.col(2) = normalise(cross(x_axis, direction_cosine.col(1)));

        trial_n = current_n = join_rows(direction_cosine, direction_cosine);

        suanpan::hacker(initial_length) = norm(x_axis);
    }

    const mat trial_disp = reshape(get_trial_displacement(element_ptr), 6, 2);
    const vec incre_disp = trial_disp.head_rows(3).col(1) - trial_disp.head_rows(3).col(0);
    const vec trial_cord = x_axis + incre_disp;

    // new length
    length = norm(trial_cord);

    elongation = dot(x_axis + trial_cord, incre_disp) / (length + initial_length); // eq. 4.98

    const mat incre_r = reshape(get_incre_displacement(element_ptr), 6, 2).eval().tail_rows(3);

    // nodal frame
    trial_n.head_cols(3) = transform::rodrigues(incre_r.col(0)) * current_n.head_cols(3);
    trial_n.tail_cols(3) = transform::rodrigues(incre_r.col(1)) * current_n.tail_cols(3);

    // reference frame
    reference = transform::rodrigues((.5 * sum(trial_disp.tail_rows(3), 1)).eval()) * direction_cosine;

    // basic deformed frame
    e(0) = normalise(trial_cord);
    for(auto I = 1u; I < 3u; ++I) e(I) = r(I) - .5 * dot(r(I), e(0)) * (r(0) + e(0));
    for(auto I = 0u; I < 3u; ++I) se(I) = transform::skew_symm(e(I));

    theta.set_size(6);
    for(auto I = 0u; I < 3u; ++I) {
        const auto J = (I + 1u) % 3u, K = (I + 2u) % 3u;
        theta(I) = asin(.5 * (dot(e(K), ni(J)) - dot(e(J), ni(K))));
        theta(I + 3llu) = asin(.5 * (dot(e(K), nj(J)) - dot(e(J), nj(K))));
    }

    const auto a = compute_a();
    const auto lr1 = compute_l(a, r(1)), lr2 = compute_l(a, r(2));
    for(auto I = 0u; I < 6u; ++I) sn(I) = transform::skew_symm(trial_n.col(I));

    const mat na = trial_n.t() * a;

    transformation.zeros(6, 12);

    // row 1
    transformation.row(0).cols(0, 2) = -e(0).t();
    transformation.row(0).cols(6, 8) = e(0).t();

    // row 2
    transformation.row(1) = ni(0).t() * lr1;
    transformation.row(1).cols(0, 2) += na.row(1);
    transformation.row(1).cols(3, 5) += e(0).t() * sni(1) - e(1).t() * sni(0);
    transformation.row(1).cols(6, 8) -= na.row(1);
    transformation.row(1) *= .5 / std::cos(theta(2));

    // row 3
    transformation.row(2) = nj(0).t() * lr1;
    transformation.row(2).cols(0, 2) += na.row(4);
    transformation.row(2).cols(6, 8) -= na.row(4);
    transformation.row(2).cols(9, 11) += e(0).t() * snj(1) - e(1).t() * snj(0);
    transformation.row(2) *= .5 / std::cos(theta(5));

    // row 4
    transformation.row(3) = -ni(0).t() * lr2;
    transformation.row(3).cols(0, 2) -= na.row(2);
    transformation.row(3).cols(3, 5) -= e(0).t() * sni(2) - e(2).t() * sni(0);
    transformation.row(3).cols(6, 8) += na.row(2);
    transformation.row(3) *= .5 / std::cos(theta(1));

    // row 5
    transformation.row(4) = -nj(0).t() * lr2;
    transformation.row(4).cols(0, 2) -= na.row(5);
    transformation.row(4).cols(6, 8) += na.row(5);
    transformation.row(4).cols(9, 11) -= e(0).t() * snj(2) - e(2).t() * snj(0);
    transformation.row(4) *= .5 / std::cos(theta(4));

    // row 6
    t6i = ni(1).t() * lr2 - ni(2).t() * lr1;
    t6j = nj(1).t() * lr2 - nj(2).t() * lr1;
    t6i.cols(3, 5) += e(1).t() * sni(2) - e(2).t() * sni(1);
    t6j.cols(9, 11) += e(1).t() * snj(2) - e(2).t() * snj(1);
    t6i *= .5 / std::cos(theta(0));
    t6j *= .5 / std::cos(theta(3));

    transformation.row(5) = t6j - t6i;
}

bool B3DC::is_nlgeom() const { return true; }

vec B3DC::to_local_vec(const vec&) const { return {elongation, theta(2), theta(5), theta(1), theta(4), theta(3) - theta(0)}; }

vec B3DC::to_global_vec(const vec& l_resistance) const { return transformation.t() * l_resistance; }

mat B3DC::to_global_geometry_mat(const mat& l_force) const {
    mat geometry(12, 12, fill::zeros);

    const auto &p1 = l_force(0), &p2 = l_force(1), &p3 = l_force(2), &p4 = l_force(3), &p5 = l_force(4), &p6 = l_force(5);

    const auto m2 = p2 * .5 / std::cos(theta(2));
    const auto m3 = p3 * .5 / std::cos(theta(5));
    const auto m4 = p4 * .5 / std::cos(theta(1));
    const auto m5 = p5 * .5 / std::cos(theta(4));
    const auto m6i = p6 * .5 / std::cos(theta(0));
    const auto m6j = p6 * .5 / std::cos(theta(3));

    const auto sa = span(0, 2), sb = span(3, 5), sc = span(6, 8), sd = span(9, 11);

    const mat a = compute_a();

    // KA
    geometry(sa, sa) = p1 * a;
    geometry(sc, sc) = geometry(sa, sa);
    geometry(sa, sc) = -geometry(sa, sa);
    geometry(sc, sa) = -geometry(sa, sa);

    // KB
    geometry += p2 * std::tan(theta(2)) * transformation.row(1).t() * transformation.row(1);
    geometry += p3 * std::tan(theta(5)) * transformation.row(2).t() * transformation.row(2);
    geometry += p4 * std::tan(theta(1)) * transformation.row(3).t() * transformation.row(3);
    geometry += p5 * std::tan(theta(4)) * transformation.row(4).t() * transformation.row(4);
    geometry -= p6 * std::tan(theta(0)) * t6i.t() * t6i;
    geometry += p6 * std::tan(theta(3)) * t6j.t() * t6j;

    // KC
    geometry += m2 * compute_g(a, r(1), ni(0)) + m3 * compute_g(a, r(1), nj(0)) - m4 * compute_g(a, r(2), ni(0)) - m5 * compute_g(a, r(2), nj(0)) + m6i * compute_g(a, r(2), nj(1)) - m6i * compute_g(a, r(1), nj(2)) - m6j * compute_g(a, r(2), ni(1)) + m6j * compute_g(a, r(1), ni(2));

    // KD
    const mat KD2 = compute_l(a, r(2)).t() * (m4 * sni(0) + m6i * sni(1)) - compute_l(a, r(1)).t() * (m2 * sni(0) + m6i * sni(2));
    const mat KD4 = compute_l(a, r(2)).t() * (m5 * snj(0) - m6j * snj(1)) - compute_l(a, r(1)).t() * (m3 * snj(0) - m6j * snj(2));

    geometry.cols(sb) += KD2;
    geometry.cols(sd) += KD4;

    // KE
    geometry.rows(sb) += KD2.t();
    geometry.rows(sd) += KD4.t();

    // KF
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

    // KF22
    geometry(sb, sb) += m2 * (se(1) * sni(0) - se(0) * sni(1)) - m4 * (se(2) * sni(0) - se(0) * sni(2)) - m6i * (se(2) * sni(1) - se(1) * sni(2));
    // KF44
    geometry(sd, sd) += m3 * (se(1) * snj(0) - se(0) * snj(1)) - m5 * (se(2) * snj(0) - se(0) * snj(2)) + m6j * (se(2) * snj(1) - se(1) * snj(2));

    return geometry;
}

mat B3DC::to_global_stiffness_mat(const mat& l_mat) const { return transformation.t() * l_mat * transformation; }
