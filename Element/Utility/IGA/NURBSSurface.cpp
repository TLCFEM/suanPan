////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2024 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "NURBSSurface.h"

vec NURBSSurface::evaluate_point(const double u, const double v, const field<vec>& polygon) const {
    const auto point = BSplineSurface::evaluate_point(u, v, polygon);

    return point.head(point.n_elem - 1) / point.back();
}

field<vec> NURBSSurface::evaluate_point_derivative(const double u, const double v, const field<vec>& polygon, sword du, sword dv) const {
    if(du < 0) du = line_u.get_order();
    if(dv < 0) dv = line_v.get_order();

    const auto point = BSplineSurface::evaluate_point_derivative(u, v, polygon, du, dv);

    const auto d = du > dv ? du : dv;
    field<vec> ders(d + 1, d + 1);

    initialize_binomial(d);

    for(auto k = 0ll; k <= du; ++k)
        for(auto l = 0ll; l <= dv; ++l) {
            ders(k, l) = point(k, l);
            for(auto i = 0ll; i <= k; ++i) for(auto j = 0ll; j <= l; ++j) if(i != 0 || j != 0) ders(k, l) -= binomial_mat(k, i) * binomial_mat(l, j) * point(i, j).back() * ders(k - i, l - j);
            ders(k, l) /= point(0, 0).back();
        }

    ders.for_each([&](vec& t_ders) { t_ders = t_ders.head(dimension - 1); });

    return ders;
}

mat NURBSSurface::evaluate_shape_function(const double u, const double v, const field<vec>& polygon) const {
    auto shape = BSplineSurface::evaluate_shape_function(u, v, polygon);

    auto sum = 0.;

    for(auto I = 0llu; I < shape.n_rows; ++I) for(auto J = 0llu; J < shape.n_cols; ++J) if(!polygon(I, J).empty()) sum += shape(I, J) * polygon(I, J).back();

    return shape / sum;
}

field<mat> NURBSSurface::evaluate_shape_function_derivative(const double u, const double v, const field<vec>& polygon, sword du, sword dv) const {
    if(du < 0) du = line_u.get_order();
    if(dv < 0) dv = line_v.get_order();

    const auto uspan = line_u.evaluate_span(u);
    const auto vspan = line_v.evaluate_span(v);

    auto uders = line_u.evaluate_basis_derivative(u, du);
    uders.resize(du + 1, uders.n_cols);
    auto vders = line_v.evaluate_basis_derivative(v, dv);
    vders.resize(dv + 1, vders.n_cols);

    const auto p = line_u.get_order();
    const auto q = line_v.get_order();

    const auto d = du > dv ? du : dv;

    initialize_binomial(d);

    field<mat> shape(d + 1, d + 1);
    shape.for_each([&](mat& t_shape) { t_shape.zeros(line_u.get_knot().n_elem - line_u.get_order() - 1, line_v.get_knot().n_elem - line_v.get_order() - 1); });

    auto sum = 0.;

    for(auto i = 0llu; i <= p; ++i) {
        const auto uind = uspan + i - p;
        for(auto j = 0llu; j <= q; ++j)
            if(const auto vind = vspan + j - q; !polygon(uind, vind).empty()) {
                sum += uders(0, i) * vders(0, j) * polygon(uind, vind).back();
                for(auto k = 0; k <= du; ++k) for(auto l = 0; l <= dv; ++l) shape(k, l)(uind, vind) = uders(k, i) * vders(l, j) * polygon(uind, vind).back();
            }
    }

    for(auto i = 0ll; i <= du; ++i)
        for(auto j = 0ll; j <= dv; ++j) {
            auto& t_shape = shape(i, j);
            for(auto k = 0ll; k <= i; ++k)
                for(auto l = 0ll; l <= j; ++l)
                    if(k != 0 || l != 0) {
                        auto weight_sum = 0.;
                        for(auto m = 0llu; m <= p; ++m) {
                            const auto uind = uspan + m - p;
                            for(auto n = 0llu; n <= q; ++n) if(const auto vind = vspan + n - q; !polygon(uind, vind).empty()) weight_sum += polygon(uind, vind).back() * uders(k, m) * vders(l, n);
                        }
                        t_shape -= binomial_mat(i, k) * binomial_mat(j, l) * weight_sum * shape(i - k, j - l);
                    }
            t_shape /= sum;
        }

    return shape;
}

NURBSSurface2D::NURBSSurface2D(vec knot_u, vec knot_v, field<vec>&& N)
    : NURBSSurface(std::move(knot_u), std::move(knot_v), 3, std::move(N)) {}

NURBSSurface3D::NURBSSurface3D(vec knot_u, vec knot_v, field<vec>&& N)
    : NURBSSurface(std::move(knot_u), std::move(knot_v), 4, std::move(N)) {}

NURBSSurface4D::NURBSSurface4D(vec knot_u, vec knot_v, field<vec>&& N)
    : NURBSSurface(std::move(knot_u), std::move(knot_v), 5, std::move(N)) {}
