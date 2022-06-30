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

#include "LinearHardeningNM.h"
#include <Toolbox/utility.h>

bool LinearHardeningNM::update_nodal_quantity(mat& jacobian, vec& residual, const double gm, const vec& q, const vec& b, const double alpha) const {
    const vec s = q - b;
    const auto h = compute_h(alpha);
    const auto f = compute_f(s, h);

    if(gm <= datum::eps && f < 0.) return false;

    jacobian.set_size(j_size, j_size);
    residual.set_size(j_size);

    const auto g = compute_df(s, h);
    const double n = norm(g);
    const vec z = g / n;
    const vec dh = -s % compute_dh(alpha) / h;

    const mat gdz = gm / n * (eye(n_size, n_size) - z * z.t()) * compute_ddf(s, h);
    const vec dgzdg = gdz * dh + z;

    residual(sa) = gm * z;
    residual(sc).fill(f);

    jacobian(sa, sa) = gdz;
    jacobian(sa, sc) = dgzdg;

    jacobian(sc, sa) = g.t();
    jacobian(sc, sc).fill(dot(g, dh));

    if(has_kinematic) {
        jacobian(sa, sb) = -gdz;
        jacobian(sc, sb) = -g.t();

        residual(sb) = -kinematic_modulus * gm * z;

        jacobian(sb, sa) = -kinematic_modulus * gdz;
        jacobian(sb, sb) = kinematic_modulus * gdz;
        jacobian(sb, sc) = -kinematic_modulus * dgzdg;
    }

    return true;
}

vec LinearHardeningNM::compute_h(const double alpha) const {
    vec h(n_size, fill::value(std::max(datum::eps, 1. + isotropic_modulus * alpha)));

    h(0) = std::max(datum::eps, 1. + 2. * isotropic_modulus * alpha);

    return h;
}

vec LinearHardeningNM::compute_dh(const double alpha) const {
    auto dh = compute_h(alpha);

    dh(0) = suanpan::approx_equal(dh(0), datum::eps) ? 0. : 2. * isotropic_modulus;

    for(auto I = 1u; I < n_size; ++I) dh(I) = suanpan::approx_equal(dh(I), datum::eps) ? 0. : isotropic_modulus;

    return dh;
}

LinearHardeningNM::LinearHardeningNM(const unsigned T, const double EEA, const double EEIS, const double HH, const double KK, const double LD, vec&& YF)
    : NonlinearNM(T, EEA, EEIS, !suanpan::approx_equal(KK, 0.), LD, std::forward<vec>(YF))
    , isotropic_modulus(HH)
    , kinematic_modulus(KK) {}

LinearHardeningNM::LinearHardeningNM(const unsigned T, const double EEA, const double EEIS, const double EEIW, const double HH, const double KK, const double LD, vec&& YF)
    : NonlinearNM(T, EEA, EEIS, EEIW, !suanpan::approx_equal(KK, 0.), LD, std::forward<vec>(YF))
    , isotropic_modulus(HH)
    , kinematic_modulus(KK) {}
