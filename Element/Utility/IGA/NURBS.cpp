////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2025 Theodore Chang
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

#include "NURBS.h"
#include "Toolbox/utility.h"

void NURBSBase::initialize_binomial(const sword d) const {
    if(static_cast<uword>(d) + 1 <= binomial_mat.n_rows) return;

    auto& t_mat = access::rw(binomial_mat);

    t_mat.set_size(d + 1, d + 1);

    for(auto i = 0ll; i <= d; ++i) for(auto j = 0ll; j <= i; ++j) t_mat(i, j) = static_cast<double>(suanpan::binomial(i, j));
}

vec NURBS::evaluate_point(const double u, const field<vec>& polygon) const {
    const auto point = BSpline::evaluate_point(u, polygon);

    return point.head(dimension - 1) / point.back();
}

field<vec> NURBS::evaluate_point_derivative(const double u, const field<vec>& polygon, sword d) const {
    if(d < 0) d = order;

    const auto point = BSpline::evaluate_point_derivative(u, polygon, d);

    initialize_binomial(d);

    field<vec> ders(d + 1);
    ders.for_each([&](vec& t_point) { t_point.zeros(dimension - 1); });

    for(auto k = 0ll; k <= d; ++k) {
        ders(k) = point(k).head(dimension - 1);
        for(auto i = 1ll; i <= k; ++i) ders(k) -= binomial_mat(k, i) * point(i).back() * ders(k - i);
        ders(k) /= point(0).back();
    }

    return ders;
}

vec NURBS::evaluate_shape_function(const double u, const field<vec>& polygon) const {
    auto shape = BSpline::evaluate_shape_function(u, polygon);

    auto sum = 0.;

    for(auto I = 0llu; I < shape.n_rows; ++I) for(auto J = 0llu; J < shape.n_cols; ++J) sum += shape(I, J) * polygon(I, J).back();

    return shape / sum;
}

field<vec> NURBS::evaluate_shape_function_derivative(const double u, const field<vec>& polygon, sword d) const {
    if(d < 0) d = order;

    const auto span = evaluate_span(u);
    auto ders = evaluate_basis_derivative(u, d);
    ders.resize(d + 1, ders.n_cols);

    initialize_binomial(d);

    field<vec> shape(d + 1);
    shape.for_each([&](vec& t_shape) { t_shape.zeros(get_number_of_control_points()); });

    auto sum = 0.;

    for(auto i = 0llu; i <= order; ++i)
        if(const auto ind = span + i - order; !polygon(ind).empty()) {
            sum += ders(0, i) * polygon(ind).back();
            for(auto j = 0; j <= d; ++j) shape(j)(ind) = ders(j, i) * polygon(ind).back();
        }

    for(auto i = 0ll; i <= d; ++i) {
        auto& t_shape = shape(i);
        for(auto j = 1ll; j <= i; ++j) {
            auto factor = 0.;
            for(auto k = 0llu; k <= order; ++k) if(const auto sind = span + k - order; !polygon(sind).empty()) factor += ders(j, k) * polygon(sind).back();
            t_shape -= binomial_mat(i, j) * factor * shape(i - j);
        }
        t_shape /= sum;
    }

    return shape;
}

NURBSCurve2D::NURBSCurve2D(vec K, field<vec>&& N)
    : NURBS(std::move(K), 3, std::move(N)) {}

NURBSCurve3D::NURBSCurve3D(vec K, field<vec>&& N)
    : NURBS(std::move(K), 4, std::move(N)) {}

NURBSCurve4D::NURBSCurve4D(vec K, field<vec>&& N)
    : NURBS(std::move(K), 5, std::move(N)) {}
