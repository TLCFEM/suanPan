////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2022 Theodore Chang
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

#include "BSplineVolume.h"

BSplineVolume::BSplineVolume(vec knot_u, vec knot_v, vec knot_w, const uword size, field<vec>&& N)
    : dimension(size)
    , net(std::forward<field<vec>>(N))
    , line_u(std::move(knot_u), size)
    , line_v(std::move(knot_v), size)
    , line_w(std::move(knot_w), size) {}

void BSplineVolume::set_control_polygon(field<vec>&& N) { net = std::forward<field<vec>>(N); }

void BSplineVolume::set_control_polygon(const field<vec>& N) { net = N; }

field<uvec> BSplineVolume::get_all_element_span() const {
    field<uvec> span(3);

    span(0) = line_u.get_all_element_span();
    span(1) = line_v.get_all_element_span();
    span(2) = line_w.get_all_element_span();

    return span;
}

uvec BSplineVolume::get_number_of_control_points() const { return {line_u.get_knot().n_elem - line_u.get_order() - 1, line_v.get_knot().n_elem - line_v.get_order() - 1, line_w.get_knot().n_elem - line_w.get_order() - 1}; }

vec BSplineVolume::evaluate_point(const double u, const double v, const double w) const { return evaluate_point(u, v, w, net); }

field<vec> BSplineVolume::evaluate_point_derivative(const double u, const double v, const double w, const sword d) const { return evaluate_point_derivative(u, v, w, net, d); }

cube BSplineVolume::evaluate_shape_function(const double u, const double v, const double w) const { return evaluate_shape_function(u, v, w, net); }

field<cube> BSplineVolume::evaluate_shape_function_derivative(const double u, const double v, const double w, const sword du, const sword dv, const sword dw) const { return evaluate_shape_function_derivative(u, v, w, net, du, dv, dw); }

vec BSplineVolume::evaluate_point(const double u, const double v, const double w, const field<vec>& polygon) const {
    const auto uspan = line_u.evaluate_span(u);
    const auto vspan = line_v.evaluate_span(v);
    const auto wspan = line_w.evaluate_span(w);

    const auto nu = line_u.evaluate_basis(u);
    const auto nv = line_v.evaluate_basis(v);
    const auto nw = line_w.evaluate_basis(w);

    const auto p = line_u.get_order();
    const auto q = line_v.get_order();
    const auto r = line_w.get_order();

    vec point(dimension, fill::zeros);

    for(auto m = 0llu; m <= r; ++m) {
        const auto wind = wspan - r + m;
        for(auto l = 0llu; l <= q; ++l) {
            const auto vind = vspan - q + l;
            for(auto k = 0llu; k <= p; ++k) if(const auto uind = uspan - p + k; !polygon(uind, vind, wind).empty()) point += nw(m) * nv(l) * nu(k) * polygon(uind, vind, wind);
        }
    }

    return point;
}

field<vec> BSplineVolume::evaluate_point_derivative(const double u, const double v, const double w, const field<vec>& polygon, const sword d) const {
    const auto p = line_u.get_order();
    const auto q = line_v.get_order();
    const auto r = line_w.get_order();

    field<vec> point(d + 1, d + 1, d + 1);

    point.for_each([&](vec& t_point) { t_point.zeros(dimension); });

    const auto du = std::min(d, static_cast<sword>(p));
    const auto dv = std::min(d, static_cast<sword>(q));
    const auto dw = std::min(d, static_cast<sword>(r));

    const auto uspan = line_u.evaluate_span(u);
    const auto vspan = line_v.evaluate_span(v);
    const auto wspan = line_w.evaluate_span(w);

    const auto nu = line_u.evaluate_basis_derivative(u, du);
    const auto nv = line_v.evaluate_basis_derivative(v, dv);
    const auto nw = line_w.evaluate_basis_derivative(w, dw);

    for(auto k = 0; k <= du; ++k)
        for(auto l = 0; l <= dv; ++l)
            for(auto m = 0; m <= dw; ++m)
                for(auto a = 0llu; a <= p; ++a) {
                    const auto uind = uspan - p + a;
                    for(auto b = 0llu; b <= q; ++b) {
                        const auto vind = vspan - q + b;
                        for(auto c = 0llu; c <= r; ++c) {
                            const auto wind = wspan - r + c;
                            if(!polygon(uind, vind, wind).empty()) point(k, l, m) += nw(m, c) * nv(l, b) * nu(k, a) * polygon(uind, vind, wind);
                        }
                    }
                }

    return point;
}

cube BSplineVolume::evaluate_shape_function(const double u, const double v, const double w, const field<vec>&) const {
    const auto uspan = line_u.evaluate_span(u);
    const auto vspan = line_v.evaluate_span(v);
    const auto wspan = line_w.evaluate_span(w);

    const auto nu = line_u.evaluate_basis(u);
    const auto nv = line_v.evaluate_basis(v);
    const auto nw = line_w.evaluate_basis(w);

    const auto p = line_u.get_order();
    const auto q = line_v.get_order();
    const auto r = line_w.get_order();

    cube point(line_u.get_knot().n_elem - p - 1, line_v.get_knot().n_elem - q - 1, line_w.get_knot().n_elem - r - 1, fill::zeros);

    for(auto m = 0llu; m <= r; ++m) {
        const auto wind = wspan - r + m;
        for(auto l = 0llu; l <= q; ++l) {
            const auto vind = vspan - q + l;
            for(auto k = 0llu; k <= p; ++k) {
                const auto uind = uspan - p + k;
                point(uind, vind, wind) = nw(m) * nv(l) * nu(k);
            }
        }
    }

    return point;
}

field<cube> BSplineVolume::evaluate_shape_function_derivative(const double u, const double v, const double w, const field<vec>&, sword du, sword dv, sword dw) const {
    const auto p = line_u.get_order();
    const auto q = line_v.get_order();
    const auto r = line_w.get_order();

    du = std::min(du, static_cast<sword>(p));
    dv = std::min(dv, static_cast<sword>(q));
    dw = std::min(dw, static_cast<sword>(r));

    const auto uspan = line_u.evaluate_span(u);
    const auto vspan = line_v.evaluate_span(v);
    const auto wspan = line_w.evaluate_span(w);

    const auto nu = line_u.evaluate_basis_derivative(u, du);
    const auto nv = line_v.evaluate_basis_derivative(v, dv);
    const auto nw = line_w.evaluate_basis_derivative(w, dw);

    field<cube> point(du + 1, dv + 1, dw + 1);

    point.for_each([&](cube& t_point) { t_point.zeros(line_u.get_knot().n_elem - p - 1, line_v.get_knot().n_elem - q - 1, line_w.get_knot().n_elem - r - 1); });

    for(auto i = 0ll; i <= du; ++i)
        for(auto j = 0ll; j <= dv; ++j)
            for(auto k = 0ll; k <= dw; ++k)
                for(auto a = 0llu; a <= p; ++a) {
                    const auto uind = uspan - p + a;
                    for(auto b = 0llu; b <= q; ++b) {
                        const auto vind = vspan - q + b;
                        for(auto c = 0llu; c <= r; ++c) {
                            const auto wind = wspan - r + c;
                            point(i, j, k)(uind, vind, wind) = nw(k, c) * nv(j, b) * nu(i, a);
                        }
                    }
                }

    return point;
}

BSplineVolume3D::BSplineVolume3D(vec knot_u, vec knot_v, vec knot_w, field<vec>&& N)
    : BSplineVolume(std::move(knot_u), std::move(knot_v), std::move(knot_w), 3, std::forward<field<vec>>(N)) {}

BSplineVolume4D::BSplineVolume4D(vec knot_u, vec knot_v, vec knot_w, field<vec>&& N)
    : BSplineVolume(std::move(knot_u), std::move(knot_v), std::move(knot_w), 4, std::forward<field<vec>>(N)) {}
