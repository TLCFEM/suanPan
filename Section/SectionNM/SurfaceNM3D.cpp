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

#include "SurfaceNM3D.h"

double SurfaceNM3D::evaluate(const double p, const double ms, const double mw, const mat& weight) const {
    double value = weight(0);

    if(weight(1) > 0.) value *= pow(p, weight(1));
    if(weight(2) > 0.) value *= pow(ms, weight(2));
    if(weight(3) > 0.) value *= pow(mw, weight(3));

    return value;
}

vec SurfaceNM3D::differentiate(const mat& weight, uword location, const uword order) {
    ++location;

    vec weight_out = weight.as_col();

    weight_out(location) = weight(location) - static_cast<double>(order);

    if(weight_out(location) < 0.) weight_out.zeros();
    else for(auto I = static_cast<uword>(weight(location)); I > static_cast<uword>(weight_out(location)); --I) weight_out(0) *= static_cast<double>(I);

    return weight_out;
}

double SurfaceNM3D::compute_sf(const vec& s, const double h) const {
    const auto p = s(0) / h;
    const auto ms = s(1) / h;
    const auto mw = s(2) / h;

    auto f = -c;
    for(auto I = 0llu; I < para_set.n_rows; ++I) f += evaluate(p, ms, mw, para_set.row(I));

    return f;
}

vec SurfaceNM3D::compute_dsf(const vec& s, const double h) const {
    const auto p = s(0) / h;
    const auto ms = s(1) / h;
    const auto mw = s(2) / h;

    vec df(3, fill::zeros);

    for(auto I = 0llu; I < para_set.n_rows; ++I) for(auto J = 0llu; J < df.n_elem; ++J) df(J) += evaluate(p, ms, mw, differentiate(para_set.row(I), J, 1));

    return df / h;
}

mat SurfaceNM3D::compute_ddsf(const vec& s, const double h) const {
    const auto p = s(0) / h;
    const auto ms = s(1) / h;
    const auto mw = s(2) / h;

    mat ddf(3, 3, fill::zeros);

    for(auto I = 0llu; I < para_set.n_rows; ++I)
        for(auto J = 0llu; J < ddf.n_rows; ++J) {
            const auto dfj = differentiate(para_set.row(I), J, 1);
            for(auto K = 0llu; K < ddf.n_cols; ++K) ddf(J, K) += evaluate(p, ms, mw, differentiate(dfj, K, 1));
        }

    return ddf * pow(h, -2.);
}

SurfaceNM3D::SurfaceNM3D(const double CC, mat&& PS)
    : para_set(PS.empty() ? mat{{1.15, 2., 0., 0.}, {1., 0., 2., 0.}, {1., 0., 0., 4.}, {3.67, 2., 2., 0.}, {3., 6., 0., 2.}, {4.65, 0., 4., 2.}} : std::forward<mat>(PS))
    , c(CC) {}
