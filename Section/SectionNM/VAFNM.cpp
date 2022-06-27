/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "VAFNM.h"
#include <Toolbox/utility.h>

bool VAFNM::update_nodal_quantity(mat& jacobian, vec& residual, const double gm, const vec& q, const vec& b, const double alpha, const vec& trial_q, const vec& bn) const {
    jacobian = eye(j_size, j_size);
    residual = zeros(j_size);

    const vec s = q - b;
    const auto h = compute_h(alpha);
    const auto f = compute_f(s, h);

    if(gm <= datum::eps && f < 0.) return false;

    const auto g = compute_df(s, h);
    const double n = norm(g);
    const vec z = g / n;
    const auto dh = -compute_dh(alpha) / h;

    const mat gdz = gm / n * (eye(n_size, n_size) - z * z.t()) * compute_ddf(s, h);
    const vec gdzdh = gdz * s * dh;

    residual(sa) = q - trial_q + gm * z;
    residual(sc).fill(f);

    jacobian(sa, sa) += gdz;
    jacobian(sa, sc) += z + gdzdh;

    jacobian(sc, sa) = g.t();
    jacobian(sc, sc).fill(dh * dot(g, s));

    if(has_kinematic) {
        jacobian(sa, sb) -= gdz;
        jacobian(sc, sb) = -g.t();

        residual(sb) = (1. + kin_base * gm) * b - bn - kin_modulus * gm * z;

        jacobian(sb, sa) -= kin_modulus * gdz;
        jacobian(sb, sb) += kin_modulus * gdz + kin_base * gm * eye(n_size, n_size);
        jacobian(sb, sc) += kin_base * b - kin_modulus * (z + gdzdh);
    }

    return true;
}

double VAFNM::compute_h(const double alpha) const { return std::max(datum::eps, 1. + iso_modulus * alpha + iso_saturation - iso_saturation * exp(-iso_decay * alpha)); }

double VAFNM::compute_dh(const double alpha) const { return suanpan::approx_equal(compute_h(alpha), datum::eps) ? 0. : iso_modulus + iso_saturation * iso_decay * exp(-iso_decay * alpha); }

VAFNM::VAFNM(const unsigned T, const double EEA, const double EEIS, const double HH, const double HS, const double HD, const double KK, const double KB, const double LD, vec&& YF)
    : NonlinearNM(T, EEA, EEIS, !suanpan::approx_equal(KK, 0.) || !suanpan::approx_equal(KB, 0.), LD, std::forward<vec>(YF))
    , iso_modulus(HH)
    , kin_modulus(KK)
    , iso_saturation(HS)
    , iso_decay(HD)
    , kin_base(KB) {}

VAFNM::VAFNM(const unsigned T, const double EEA, const double EEIS, const double EEIW, const double HH, const double HS, const double HD, const double KK, const double KB, const double LD, vec&& YF)
    : NonlinearNM(T, EEA, EEIS, EEIW, !suanpan::approx_equal(KK, 0.) || !suanpan::approx_equal(KB, 0.), LD, std::forward<vec>(YF))
    , iso_modulus(HH)
    , kin_modulus(KK)
    , iso_saturation(HS)
    , iso_decay(HD)
    , kin_base(KB) {}
