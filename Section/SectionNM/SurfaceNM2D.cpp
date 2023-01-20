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

#include "SurfaceNM2D.h"

double SurfaceNM2D::evaluate(const double p, const double ms, const mat& weight) {
    double value = weight(0);

    if(static_cast<int>(weight(1)) > 0) value *= pow(p, weight(1));
    if(static_cast<int>(weight(2)) > 0) value *= pow(ms, weight(2));

    return value;
}

vec SurfaceNM2D::differentiate(const mat& weight, uword location, const uword order) {
    ++location;

    vec weight_out = weight.as_col();

    weight_out(location) = weight(location) - static_cast<double>(order);

    if(static_cast<int>(weight_out(location)) < 0) weight_out.zeros();
    else for(auto I = static_cast<uword>(weight(location)); I > static_cast<uword>(weight_out(location)); --I) weight_out(0) *= static_cast<double>(I);

    return weight_out;
}

double SurfaceNM2D::compute_sf(const vec& s, const vec& h) const {
    const auto p = s(0) / h(0);
    const auto ms = s(1) / h(1);

    auto f = -c;
    for(auto I = 0llu; I < para_set.n_rows; ++I) f += evaluate(p, ms, para_set.row(I));

    return f;
}

vec SurfaceNM2D::compute_dsf(const vec& s, const vec& h) const {
    const auto p = s(0) / h(0);
    const auto ms = s(1) / h(1);

    vec df(2, fill::zeros);

    for(auto I = 0llu; I < para_set.n_rows; ++I) for(auto J = 0llu; J < df.n_elem; ++J) df(J) += evaluate(p, ms, differentiate(para_set.row(I), J, 1));

    return df / h;
}

mat SurfaceNM2D::compute_ddsf(const vec& s, const vec& h) const {
    const auto p = s(0) / h(0);
    const auto ms = s(1) / h(1);

    mat ddf(2, 2, fill::zeros);

    for(auto I = 0llu; I < para_set.n_rows; ++I)
        for(auto J = 0llu; J < ddf.n_rows; ++J) {
            const auto dfj = differentiate(para_set.row(I), J, 1);
            for(auto K = 0llu; K < ddf.n_cols; ++K) ddf(J, K) += evaluate(p, ms, differentiate(dfj, K, 1));
        }

    return ddf / (h * h.t());
}

SurfaceNM2D::SurfaceNM2D(const double CC, mat&& PS)
    : para_set(PS.empty() ? mat{{1.15, 2., 0.}, {1., 0., 2.}, {3.67, 2., 2.}} : std::forward<mat>(PS))
    , c(CC) {}
