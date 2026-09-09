/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

pod6 CDP::compute_backbone(const double kappa, const double a_tc, const double cb_tc, const double f_tc) {
    static constexpr auto SLOPE_BOUND = 1e2;

    pod6 out;

    const auto s_phi = std::sqrt(1. + a_tc * (a_tc + 2.) * kappa);
    const auto t_phi = (1. + .5 * a_tc) / s_phi;
    const auto b_phi = (1. + a_tc - s_phi) / a_tc;

    out[1] = f_tc * s_phi * b_phi;                    // f
    out[4] = f_tc * t_phi * (1. + a_tc - 2. * s_phi); // \md{f}

    if(const auto p_phi = std::pow(b_phi, cb_tc), segment = SLOPE_BOUND - SLOPE_BOUND * kappa; p_phi > segment) {
        out[0] = 1. - segment;
        out[3] = SLOPE_BOUND;
    }
    else {
        out[0] = 1. - p_phi;                    // original d
        out[3] = cb_tc * t_phi * p_phi / b_phi; // \md{d}
    }

    out[2] = out[1] / (1. - out[0]);                     // \bar{f}
    out[5] = (out[4] + out[2] * out[3]) / (1. - out[0]); // \md{\bar{f}}

    return out;
}

pod6 CDP::compute_tension_backbone(const double kappa) const { return compute_backbone(kappa, a_t, cb_t, f_t); }

pod6 CDP::compute_compression_backbone(const double kappa) const { return compute_backbone(kappa, a_c, cb_c, f_c); }

CDP::CDP(const bool CHECK_INPUT, const unsigned T, const double E, const double V, const double ST, const double SC, const double GT, const double GC, const double AT, const double AC, const double DT, const double DC, const double AP, const double BC, const double S, const double R)
    : NonlinearCDP(T, E, V, GT, GC, AP, BC, S, R)
    , a_t(AT < 1. ? AT : .5)
    , cb_t(0.)
    , f_t(std::fabs(ST))
    , a_c(AC > 1. ? AC : 4.)
    , cb_c(0.)
    , f_c(-std::fabs(SC) * 4. * a_c * std::pow(1. + a_c, -2.)) {
    // tension
    double ratio_t{};
    if(CHECK_INPUT) {
        const auto half_elastic_strain = .5 * f_t / elastic_modulus;
        // try to find the total strain at half stress point, eq. (20), eq. (24), eq. (27) and eq. (49)
        // this defines a minimum stiffness degradation
        ratio_t = half_elastic_strain / (std::log(1. + a_t + std::sqrt(1. + a_t * a_t)) / f_t / (1. + .5 * a_t) * g_t + half_elastic_strain);
        if(ratio_t >= DT) suanpan_warning("A minimum tension degradation of {:.3f} is required, resetting.\n", ratio_t);
    }
    access::rw(cb_t) = std::log(std::clamp(DT, ratio_t, 1. - datum::eps)) / std::log(.5 * (1. + a_t - std::sqrt(1. + a_t * a_t)) / a_t);

    // compression
    double ratio_c{};
    if(CHECK_INPUT) {
        const auto peak_elastic_strain = -std::fabs(SC) / elastic_modulus;
        // try to find the total strain at peak stress point, eq. (20), eq. (24), eq. (27) and eq. (46)
        // this defines a minimum stiffness degradation
        ratio_c = peak_elastic_strain / (std::log(2. * a_c / (1. + a_c)) / f_c / (1. + .5 * a_c) * g_c + peak_elastic_strain);
        if(ratio_c >= DC) suanpan_warning("A minimum compression degradation of {:.3f} is required, resetting.\n", ratio_c);
    }
    access::rw(cb_c) = std::log(std::clamp(DC, ratio_c, 1. - datum::eps)) / std::log(.5 + .5 / a_c);

    if(a_t < 1.) {
        const auto min_g_t = f_t * f_t * (1. + .5 * a_t) * (1. - a_t) / elastic_modulus;
        if(g_t < min_g_t) suanpan_warning("There is a risk of snap-back, the tension energy must be greater than {:.4e}.\n", min_g_t);
    }
    if(a_c < 1.) {
        const auto min_g_c = f_c * f_c * (1. + .5 * a_c) * (1. - a_c) / elastic_modulus;
        if(g_c < min_g_c) suanpan_warning("There is a risk of snap-back, the compression energy must be greater than {:.4e}.\n", min_g_c);
    }
}

unique_ptr<Material> CDP::unique_copy() { return std::make_unique<CDP>(*this); }
