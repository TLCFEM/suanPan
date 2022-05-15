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

#include "CDP.h"

podarray<double> CDP::compute_tension_backbone(const double kappa) const {
    podarray<double> out(6);

    const auto s_phi = sqrt(1. + a_t * (a_t + 2.) * kappa);
    const auto t_phi = (1. + .5 * a_t) / s_phi;
    const auto b_phi = (1. + a_t - s_phi) / a_t;
    const auto p_phi = pow(b_phi, cb_t);

    out(0) = 1. - p_phi;                                    // d
    out(1) = f_t * s_phi * b_phi;                           // f
    out(2) = out(1) / p_phi;                                // \bar{f}
    out(3) = cb_t * t_phi * p_phi / b_phi;                  // \md{d}
    out(4) = f_t * t_phi * (1. + a_t - 2. * s_phi);         // \md{f}
    out(5) = (out(4) + f_t * t_phi * cb_t * s_phi) / p_phi; // \md{\bar{f}}

    return out;
}

podarray<double> CDP::compute_compression_backbone(const double kappa) const {
    podarray<double> out(6);

    const auto s_phi = sqrt(1. + a_c * (a_c + 2.) * kappa);
    const auto t_phi = (1. + .5 * a_c) / s_phi;
    const auto b_phi = (1. + a_c - s_phi) / a_c;
    const auto p_phi = pow(b_phi, cb_c);

    out(0) = 1. - p_phi;                                    // d
    out(1) = f_c * s_phi * b_phi;                           // f
    out(2) = out(1) / p_phi;                                // \bar{f}
    out(3) = cb_c * t_phi * p_phi / b_phi;                  // \md{d}
    out(4) = f_c * t_phi * (1. + a_c - 2. * s_phi);         // \md{f}
    out(5) = (out(4) + f_c * t_phi * cb_c * s_phi) / p_phi; // \md{\bar{f}}

    return out;
}

CDP::CDP(const unsigned T, const double E, const double V, const double ST, const double SC, const double GT, const double GC, const double AT, const double AC, const double DT, const double DC, const double AP, const double BC, const double S, const double R)
    : NonlinearCDP(T, E, V, GT, GC, AP, BC, S, R)
    , a_t(AT < 1. ? AT : .5)
    , cb_t(log(DT < 1. ? DT : .95) / log(.5 * (1. + a_t - sqrt(1. + a_t * a_t)) / a_t))
    , f_t(fabs(ST))
    , a_c(AC > 1. ? AC : 4.)
    , cb_c(log(DC < 1. ? DC : .95) / log(.5 + .5 / a_c))
    , f_c(-fabs(SC) * 4. * a_c * pow(1. + a_c, -2.)) {
    // tension
    const auto half_stress = .5 * f_t;
    const auto half_strain = log(1. + a_t + sqrt(1. + a_t * a_t)) / f_t / (1. + .5 * a_t) * g_t + half_stress / elastic_modulus;
    if(const auto ratio_t = half_stress / half_strain / elastic_modulus; ratio_t >= DT) {
        suanpan_warning("CDP: requires a minimum tension degradation of %.2f, now reset it.\n", ratio_t);
        access::rw(cb_t) = log(ratio_t >= .9 ? ratio_t : ratio_t + .05) / log(.5 * (1. + a_t - sqrt(1. + a_t * a_t)) / a_t);
    }
    // compression
    const auto peak_stress = -fabs(SC);
    const auto peak_strain = log(2. * a_c / (1. + a_c)) / f_c / (1. + .5 * a_c) * g_c + peak_stress / elastic_modulus;
    if(const auto ratio_c = peak_stress / peak_strain / elastic_modulus; ratio_c >= DC) {
        suanpan_warning("CDP: requires a minimum compression degradation of %.2f, now reset it.\n", ratio_c);
        access::rw(cb_c) = log(ratio_c >= .9 ? ratio_c : ratio_c + .05) / log(.5 + .5 / a_c);
    }
}

unique_ptr<Material> CDP::get_copy() { return make_unique<CDP>(*this); }
