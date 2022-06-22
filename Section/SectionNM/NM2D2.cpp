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

// ReSharper disable IdentifierTypo
#include "NM2D2.h"

double NM2D2::evaluate(const double p, const double ms, const mat& weight) const {
    double value = weight(0);

    if(weight(1) > 0.) value *= pow(p, weight(1));
    if(weight(2) > 0.) value *= pow(ms, weight(2));

    return value;
}

vec NM2D2::differentiate(const mat& weight, uword location, const uword order) {
    ++location;

    vec weight_out = weight.as_col();

    weight_out(location) = weight(location) - static_cast<double>(order);

    if(weight_out(location) < 0.) weight_out.zeros();
    else for(auto I = static_cast<uword>(weight(location)); I > static_cast<uword>(weight_out(location)); --I) weight_out(0) *= static_cast<double>(I);

    return weight_out;
}

double NM2D2::compute_f(const vec& s, const double alpha) const {
    const vec iso_force = yield_force * std::max(datum::eps, 1. + h * alpha);

    const auto p = s(0) / iso_force(0);
    const auto ms = s(1) / iso_force(1);

    auto f = -c;
    for(auto I = 0llu; I < para_set.n_rows; ++I) f += evaluate(p, ms, para_set.row(I));

    return f;
}

double NM2D2::compute_dh(const vec& s, const double alpha) const {
    const auto iso_factor = std::max(datum::eps, 1. + h * alpha);
    const vec iso_force = yield_force * iso_factor;

    const auto p = s(0) / iso_force(0);
    const auto ms = s(1) / iso_force(1);

    vec df(2, fill::zeros);

    for(auto I = 0llu; I < para_set.n_rows; ++I) for(auto J = 0llu; J < df.n_elem; ++J) df(J) += evaluate(p, ms, differentiate(para_set.row(I), J, 1));

    return dot(df, -h * pow(iso_factor, -2.) * s / yield_force);
}

vec NM2D2::compute_df(const vec& s, const double alpha) const {
    const vec iso_force = yield_force * std::max(datum::eps, 1. + h * alpha);

    const auto p = s(0) / iso_force(0);
    const auto ms = s(1) / iso_force(1);

    vec df(2, fill::zeros);

    for(auto I = 0llu; I < para_set.n_rows; ++I) for(auto J = 0llu; J < df.n_elem; ++J) df(J) += evaluate(p, ms, differentiate(para_set.row(I), J, 1));

    return df / iso_force;
}

mat NM2D2::compute_ddf(const vec& s, const double alpha) const {
    const vec iso_force = yield_force * std::max(datum::eps, 1. + h * alpha);

    const auto p = s(0) / iso_force(0);
    const auto ms = s(1) / iso_force(1);

    mat ddf(2, 2, fill::zeros);

    for(auto I = 0llu; I < para_set.n_rows; ++I)
        for(auto J = 0llu; J < ddf.n_rows; ++J) {
            const auto dfj = differentiate(para_set.row(I), J, 1);
            for(auto K = 0llu; K < ddf.n_cols; ++K) ddf(J, K) += evaluate(p, ms, differentiate(dfj, K, 1));
        }

    return ddf / (iso_force * iso_force.t());
}

NM2D2::NM2D2(const unsigned T, const double EEA, const double EEIS, const double NP, const double MSP, const double CC, const double HH, const double KK, const double LD, mat&& PS)
    : NonlinearNM(T, EEA, EEIS, KK, LD)
    , para_set(PS.empty() ? mat{{1.15, 2., 0.}, {1., 0., 2.}, {3.67, 2., 2.}} : std::forward<mat>(PS))
    , yield_force{NP, MSP}
    , c(CC)
    , h(HH) {}

unique_ptr<Section> NM2D2::get_copy() { return make_unique<NM2D2>(*this); }

void NM2D2::print() {
    suanpan_info("A N-M based section. doi: 10.1007/978-94-007-6573-3_3\n");
    current_deformation.t().print("Local Deformation:");
    current_resistance.t().print("Local Resistance:");
}
