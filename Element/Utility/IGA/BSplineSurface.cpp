////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2023 Theodore Chang
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

#include "BSplineSurface.h"

BSplineSurface::BSplineSurface(vec knot_u, vec knot_v, const uword size, field<vec>&& N)
    : dimension(size)
    , net(std::move(N))
    , line_u(std::move(knot_u), size)
    , line_v(std::move(knot_v), size) {}

void BSplineSurface::set_control_polygon(field<vec>&& N) { net = std::move(N); }

void BSplineSurface::set_control_polygon(const field<vec>& N) { net = N; }

field<uvec> BSplineSurface::get_all_element_span() const {
    field<uvec> span(2);

    span(0) = line_u.get_all_element_span();
    span(1) = line_v.get_all_element_span();

    return span;
}

uvec BSplineSurface::get_number_of_control_points() const { return {line_u.get_knot().n_elem - line_u.get_order() - 1, line_v.get_knot().n_elem - line_v.get_order() - 1}; }

vec BSplineSurface::evaluate_point(const double u, const double v) const { return evaluate_point(u, v, net); }

field<vec> BSplineSurface::evaluate_point_derivative(const double u, const double v, const sword du, const sword dv) const { return evaluate_point_derivative(u, v, net, du, dv); }

mat BSplineSurface::evaluate_shape_function(const double u, const double v) const { return evaluate_shape_function(u, v, net); }

field<mat> BSplineSurface::evaluate_shape_function_derivative(const double u, const double v, const sword du, const sword dv) const { return evaluate_shape_function_derivative(u, v, net, du, dv); }

vec BSplineSurface::evaluate_point(const double u, const double v, const field<vec>& polygon) const {
    const auto uspan = line_u.evaluate_span(u);
    const auto vspan = line_v.evaluate_span(v);

    const auto nu = line_u.evaluate_basis(u);
    const auto nv = line_v.evaluate_basis(v);

    const auto p = line_u.get_order();
    const auto q = line_v.get_order();

    vec point(dimension, fill::zeros);

    for(auto k = 0llu; k <= p; ++k) {
        const auto uind = uspan - p + k;
        for(auto l = 0llu; l <= q; ++l) if(const auto vind = vspan - q + l; !polygon(uind, vind).empty()) point += nu(k) * nv(l) * polygon(uind, vind);
    }

    return point;
}

field<vec> BSplineSurface::evaluate_point_derivative(const double u, const double v, const field<vec>& polygon, sword du, sword dv) const {
    const auto p = line_u.get_order();
    const auto q = line_v.get_order();

    if(du < 0) du = p;
    if(dv < 0) dv = q;

    field<vec> point(du + 1, dv + 1);

    point.for_each([&](vec& t_point) { t_point.zeros(dimension); });

    du = std::min(du, static_cast<sword>(p));
    dv = std::min(dv, static_cast<sword>(q));

    const auto uspan = line_u.evaluate_span(u);
    const auto vspan = line_v.evaluate_span(v);

    const auto nu = line_u.evaluate_basis_derivative(u, du);
    const auto nv = line_v.evaluate_basis_derivative(v, dv);

    for(auto k = 0; k <= du; ++k)
        for(auto l = 0; l <= dv; ++l)
            for(auto r = 0llu; r <= p; ++r) {
                const auto uind = uspan - p + r;
                for(auto s = 0llu; s <= q; ++s) if(const auto vind = vspan - q + s; !polygon(uind, vind).empty()) point(k, l) += nu(k, r) * nv(l, s) * polygon(uind, vind);
            }

    return point;
}

mat BSplineSurface::evaluate_shape_function(const double u, const double v, const field<vec>&) const {
    const auto uspan = line_u.evaluate_span(u);
    const auto vspan = line_v.evaluate_span(v);

    const auto nu = line_u.evaluate_basis(u);
    const auto nv = line_v.evaluate_basis(v);

    const auto p = line_u.get_order();
    const auto q = line_v.get_order();

    mat shape(line_u.get_knot().n_elem - p - 1, line_v.get_knot().n_elem - q - 1, fill::zeros);

    for(auto k = 0llu; k <= p; ++k) {
        const auto uind = uspan - p + k;
        for(auto l = 0llu; l <= q; ++l) {
            const auto vind = vspan - q + l;
            shape(uind, vind) = nu(k) * nv(l);
        }
    }

    return shape;
}

field<mat> BSplineSurface::evaluate_shape_function_derivative(const double u, const double v, const field<vec>&, sword du, sword dv) const {
    const auto p = line_u.get_order();
    const auto q = line_v.get_order();

    const auto uspan = line_u.evaluate_span(u);
    const auto vspan = line_v.evaluate_span(v);

    du = du < 0 ? p : std::min(du, static_cast<sword>(p));
    dv = dv < 0 ? q : std::min(dv, static_cast<sword>(q));

    const auto nu = line_u.evaluate_basis_derivative(u, du);
    const auto nv = line_v.evaluate_basis_derivative(v, dv);

    field<mat> shape(du + 1, dv + 1);
    shape.for_each([&](mat& t_shape) { t_shape.zeros(line_u.get_knot().n_elem - p - 1, line_v.get_knot().n_elem - q - 1); });

    for(auto i = 0ll; i <= du; ++i)
        for(auto j = 0ll; j <= dv; ++j)
            for(auto k = 0llu; k <= p; ++k) {
                const auto uind = uspan - p + k;
                for(auto l = 0llu; l <= q; ++l) {
                    const auto vind = vspan - q + l;
                    shape(i, j)(uind, vind) = nu(i, k) * nv(j, l);
                }
            }

    return shape;
}

BSplineSurface2D::BSplineSurface2D(vec knot_u, vec knot_v, field<vec>&& N)
    : BSplineSurface(std::move(knot_u), std::move(knot_v), 2, std::move(N)) {}

BSplineSurface3D::BSplineSurface3D(vec knot_u, vec knot_v, field<vec>&& N)
    : BSplineSurface(std::move(knot_u), std::move(knot_v), 3, std::move(N)) {}

BSplineSurface4D::BSplineSurface4D(vec knot_u, vec knot_v, field<vec>&& N)
    : BSplineSurface(std::move(knot_u), std::move(knot_v), 4, std::move(N)) {}
