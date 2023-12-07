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

#include "NURBSVolume.h"

vec NURBSVolume::evaluate_point(const double u, const double v, const double w, const field<vec>& polygon) const {
    const auto point = BSplineVolume::evaluate_point(u, v, w, polygon);

    return point.head(point.n_elem - 1) / point.back();
}

field<vec> NURBSVolume::evaluate_point_derivative(const double u, const double v, const double w, const field<vec>& polygon, const sword d) const {
    const auto point = BSplineVolume::evaluate_point_derivative(u, v, w, polygon, d);

    field<vec> ders(d + 1, d + 1, d + 1);
    // ders.for_each([&](vec& t_ders) { t_ders.zeros(dimension); });

    initialize_binomial(d);

    for(auto l = 0ll; l <= d; ++l) {
        for(auto m = 0ll; m <= d; ++m) {
            for(auto n = 0ll; n <= d; ++n) {
                auto& t_ders = ders(l, m, n) = point(l, m, n);
                for(auto i = 0ll; i <= l; ++i) for(auto j = 0ll; j <= m; ++j) for(auto k = 0ll; k <= n; ++k) if(i != 0 || j != 0 || k != 0) t_ders -= binomial_mat(l, i) * binomial_mat(m, j) * binomial_mat(n, k) * point(i, j, k).back() * ders(l - i, m - j, n - k);
                t_ders /= point(0, 0, 0).back();
            }
        }
    }

    ders.for_each([&](vec& t_ders) { t_ders = t_ders.head(dimension - 1); });

    return ders;
}

cube NURBSVolume::evaluate_shape_function(const double u, const double v, const double w, const field<vec>& polygon) const {
    auto shape = BSplineVolume::evaluate_shape_function(u, v, w, polygon);

    auto sum = 0.;

    for(auto I = 0llu; I < shape.n_rows; ++I) for(auto J = 0llu; J < shape.n_cols; ++J) for(auto K = 0llu; K < shape.n_slices; ++K) if(!polygon(I, J, K).empty()) sum += shape(I, J, K) * polygon(I, J, K).back();

    return shape / sum;
}

field<cube> NURBSVolume::evaluate_shape_function_derivative(const double u, const double v, const double w, const field<vec>& polygon, sword du, sword dv, sword dw) const {
    if(du < 0) du = line_u.get_order();
    if(dv < 0) dv = line_v.get_order();
    if(dw < 0) dw = line_w.get_order();

    const auto uspan = line_u.evaluate_span(u);
    const auto vspan = line_v.evaluate_span(v);
    const auto wspan = line_w.evaluate_span(w);

    auto uders = line_u.evaluate_basis_derivative(u, du);
    uders.resize(du + 1, uders.n_cols);
    auto vders = line_v.evaluate_basis_derivative(v, dv);
    vders.resize(dv + 1, vders.n_cols);
    auto wders = line_w.evaluate_basis_derivative(w, dw);
    wders.resize(dw + 1, wders.n_cols);

    const auto p = line_u.get_order();
    const auto q = line_v.get_order();
    const auto r = line_w.get_order();

    auto d = du;
    if(dv > d) d = dv;
    if(dw > d) d = dw;

    initialize_binomial(d);

    field<cube> shape(d + 1, d + 1, d + 1);
    shape.for_each([&](cube& t_shape) { t_shape.zeros(line_u.get_knot().n_elem - line_u.get_order() - 1, line_v.get_knot().n_elem - line_v.get_order() - 1, line_w.get_knot().n_elem - line_w.get_order() - 1); });

    auto sum = 0.;

    for(auto i = 0llu; i <= p; ++i) {
        const auto uind = uspan + i - p;
        for(auto j = 0llu; j <= q; ++j) {
            const auto vind = vspan + j - q;
            for(auto k = 0llu; k <= r; ++k)
                if(const auto wind = wspan + k - r; !polygon(uind, vind, wind).empty()) {
                    sum += uders(0, i) * vders(0, j) * wders(0, k) * polygon(uind, vind, wind).back();
                    for(auto l = 0; l <= du; ++l) for(auto m = 0; m <= dv; ++m) for(auto n = 0; n <= dw; ++n) shape(l, m, n)(uind, vind, wind) = uders(l, i) * vders(m, j) * wders(n, k) * polygon(uind, vind, wind).back();
                }
        }
    }

    for(auto i = 0ll; i <= du; ++i)
        for(auto j = 0ll; j <= dv; ++j)
            for(auto k = 0ll; k <= dw; ++k) {
                auto& t_shape = shape(i, j, k);
                for(auto l = 0ll; l <= i; ++l)
                    for(auto m = 0ll; m <= j; ++m)
                        for(auto n = 0ll; n <= k; ++n)
                            if(l != 0 || m != 0 || n != 0) {
                                auto weight_sum = 0.;
                                for(auto x = 0llu; x <= p; ++x) {
                                    const auto uind = uspan + x - p;
                                    for(auto y = 0llu; y <= q; ++y) {
                                        const auto vind = vspan + y - q;
                                        for(auto z = 0llu; z <= r; ++z) {
                                            const auto wind = wspan + z - r;
                                            if(!polygon(uind, vind, wind).empty()) weight_sum += polygon(uind, vind, wind).back() * uders(l, x) * vders(m, y) * wders(n, z);
                                        }
                                    }
                                }
                                t_shape -= binomial_mat(i, l) * binomial_mat(j, m) * binomial_mat(k, n) * weight_sum * shape(i - l, j - m, k - n);
                            }
                t_shape /= sum;
            }

    return shape;
}

NURBSVolume3D::NURBSVolume3D(vec knot_u, vec knot_v, vec knot_w, field<vec>&& N)
    : NURBSVolume(std::move(knot_u), std::move(knot_v), std::move(knot_w), 4, std::move(N)) {}

NURBSVolume4D::NURBSVolume4D(vec knot_u, vec knot_v, vec knot_w, field<vec>&& N)
    : NURBSVolume(std::move(knot_u), std::move(knot_v), std::move(knot_w), 5, std::move(N)) {}
