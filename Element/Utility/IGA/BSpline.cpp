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

#include "BSpline.h"

void IGA::convert_to_weighted(mat& polygon) { for(auto I = 0llu; I < polygon.n_rows - 1; ++I) polygon.row(I) %= polygon.tail_rows(1); }

void IGA::convert_to_weighted(field<vec>& polygon) {
    polygon.for_each([](vec& t_point) {
        if(t_point.empty()) return;
        t_point.head(t_point.n_elem - 1) *= t_point.back();
    });
}

uword IGA::compute_order(const vec& knot) {
    const auto b_order = find(knot.min() == knot).eval().n_elem - 1;

    if(const auto e_order = find(knot.max() == knot).eval().n_elem - 1; b_order == e_order) return b_order;

    throw std::invalid_argument("inconsistent order detected");
}

uword IGA::compute_number_of_elements(const vec& knot) { return unique(knot).eval().n_elem - 1; }

uword IGA::compute_number_of_control_points(const vec& knot) { return knot.n_elem - compute_order(knot) - 1; }

uvec IGA::compute_all_element_span(const vec& knot) {
    std::vector<uword> span;
    span.reserve(knot.n_elem);

    for(auto I = 0llu, J = 1llu; I < knot.n_elem - 1; ++I, ++J) if(knot(I) < knot(J)) span.emplace_back(I);

    return span;
}

BSpline::BSpline(vec K, const uword S, field<vec>&& N)
    : dimension(S)
    , knot(std::move(K))
    , net(std::move(N)) {}

uword BSpline::get_order() const { return order; }

uword BSpline::get_number_of_control_points() const { return knot.n_elem - order - 1; }

void BSpline::set_control_polygon(field<vec>&& N) { net = std::move(N); }

void BSpline::set_control_polygon(const field<vec>& N) { net = N; }

const vec& BSpline::get_knot() const { return knot; }

uvec BSpline::get_all_element_span() const { return IGA::compute_all_element_span(knot); }

/**
 * \brief Algorithm A2.1
 * \param u parameter
 * \return span index
 */
uword BSpline::evaluate_span(const double u) const {
    const auto n = knot.n_elem - order - 2llu;

    if(fabs(u - knot(n + 1)) <= 1E-14) return n;

    auto low = order, high = n + 1;

    uword mid;
    while(true) {
        mid = (low + high) / 2;
        if(u < knot(mid)) high = mid;
        else if(u >= knot(mid + 1)) low = mid;
        else break;
    }

    return mid;
}

/**
 * \brief Algorithm A2.2
 * \param u parameter
 * \param p degree
 * \return basis function
 */
vec BSpline::evaluate_basis(const double u, sword p) const {
    p = p < 0 ? static_cast<sword>(order) : std::min(p, static_cast<sword>(order));

    ++p;

    vec basis(p, fill::ones), left(p), right(p);

    const auto i = evaluate_span(u);

    for(auto j = 1ll; j < p; ++j) {
        left(j) = u - knot(i + 1 - j);
        right(j) = knot(i + j) - u;
        basis(j) = 0.;
        for(auto r = 0ll; r < j; ++r) {
            const auto &c_right = right(r + 1), &c_left = left(j - r);
            const auto factor = basis(r) / (c_right + c_left);
            basis(r) = basis(j) + c_right * factor;
            basis(j) = c_left * factor;
        }
    }

    return basis;
}

/**
 * \brief Algorithm A2.3
 * \param u parameter
 * \param n order
 * \param p degree
 * \return derivative of basis function
 */
mat BSpline::evaluate_basis_derivative(const double u, sword n, sword p) const {
    p = p < 0 ? static_cast<sword>(order) : std::min(p, static_cast<sword>(order));
    n = n < 0 ? 0 : std::min(n, p);

    const auto i = evaluate_span(u);

    const auto pp = p + 1;

    mat ndu(pp, pp), ders(n + 1, pp, fill::zeros), a(2, pp);

    vec left(pp), right(pp);

    ndu(0, 0) = 1.;
    for(auto j = 1ll; j <= p; ++j) {
        left(j) = u - knot(i + 1 - j);
        right(j) = knot(i + j) - u;
        ndu(j, j) = 0.;
        for(auto r = 0ll; r < j; ++r) {
            const auto &c_right = right(r + 1), &c_left = left(j - r);
            ndu(j, r) = c_right + c_left;
            const auto factor = ndu(r, j - 1) / ndu(j, r);
            ndu(r, j) = ndu(j, j) + c_right * factor;
            ndu(j, j) = c_left * factor;
        }
    }

    for(auto j = 0; j <= p; ++j) ders(0, j) = ndu(j, p);

    for(auto r = 0ll; r <= p; ++r) {
        // do not optimise due to swap
        // ReSharper disable once CppTooWideScope
        auto s1 = 0;
        // ReSharper disable once CppTooWideScope
        auto s2 = 1;
        a(0, 0) = 1.;
        for(auto k = 1ll; k <= n; ++k) {
            const auto rk = r - k, pk = p - k;
            auto& d = ders(k, r);
            if(r >= k) d += (a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk)) * ndu(rk, pk);
            for(auto j = rk >= -1 ? 1 : -rk; j <= (r <= pk + 1 ? k - 1 : p - r); ++j) d += (a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j)) * ndu(rk + j, pk);
            if(r <= pk) d += (a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r)) * ndu(r, pk);
            std::swap(s1, s2);
        }
    }

    auto r = static_cast<double>(p);
    for(auto k = 1; k <= n; ++k) {
        for(auto j = 0; j <= p; ++j) ders(k, j) *= r;
        r *= static_cast<double>(p - k);
    }

    return ders;
}

vec BSpline::evaluate_point(const double u) const { return evaluate_point(u, net); }

field<vec> BSpline::evaluate_point_derivative(const double u, const sword du) const { return evaluate_point_derivative(u, net, du); }

vec BSpline::evaluate_shape_function(const double u) const { return evaluate_shape_function(u, net); }

field<vec> BSpline::evaluate_shape_function_derivative(const double u, const sword du) const { return evaluate_shape_function_derivative(u, net, du); }

/**
 * \brief Algorithm A3.1
 * \param u parameter
 * \param polygon control net polygon
 * \return vector contains location of point
 */
vec BSpline::evaluate_point(const double u, const field<vec>& polygon) const {
    const auto span = evaluate_span(u);
    const auto basis = evaluate_basis(u);

    suanpan_assert([&] { if(polygon.n_elem < knot.n_elem - order - 1) throw invalid_argument("need more control points"); });

    vec point(dimension, fill::zeros);

    for(auto i = 0llu; i <= order; ++i) if(const auto sind = span + i - order; !polygon(sind).empty()) point += basis(i) * polygon(sind);

    return point;
}

/**
 * \brief Algorithm A3.2
 * \param u parameter
 * \param polygon control net polygon
 * \param du degree
 * \return matrix contains derivatives arranged in the following layout
 *         | 0th derivative | 1st derivative | 2nd derivative | ... |
 *         | x component    | x component    | x component    | ... |
 *         | y component    | y component    | y component    | ... |
 *         | ...            | ...            | ...            | ... |
 *
 */
field<vec> BSpline::evaluate_point_derivative(const double u, const field<vec>& polygon, sword du) const {
    if(du < 0) du = order;

    field<vec> point(du + 1);
    point.for_each([&](vec& t_point) { t_point.zeros(dimension); });

    du = std::min(du, static_cast<sword>(order));

    const auto span = evaluate_span(u);
    const auto ders = evaluate_basis_derivative(u, du);

    for(auto k = 0ll; k <= du; ++k) for(auto j = 0llu; j <= order; ++j) if(const auto sind = span + j - order; !polygon(sind).empty()) point(k) += ders(k, j) * polygon(sind);

    return point;
}

vec BSpline::evaluate_shape_function(const double u, const field<vec>&) const {
    vec shape(get_number_of_control_points(), fill::zeros);

    const auto span = evaluate_span(u);
    const auto basis = evaluate_basis(u);

    for(auto i = 0llu; i <= order; ++i) shape(span + i - order) = basis(i);

    return shape;
}

field<vec> BSpline::evaluate_shape_function_derivative(const double u, const field<vec>&, sword du) const {
    du = du < 0 ? order : std::min(du, static_cast<sword>(order));

    field<vec> shape(du + 1);
    shape.for_each([&](vec& t_shape) { t_shape.zeros(get_number_of_control_points()); });

    const auto span = evaluate_span(u);
    const auto ders = evaluate_basis_derivative(u, du);

    for(auto i = 0ll; i <= du; ++i) for(auto j = 0llu; j <= order; ++j) shape(i)(span + j - order) = ders(i, j);

    return shape;
}

BSplineCurve2D::BSplineCurve2D(vec K, field<vec>&& N)
    : BSpline(std::move(K), 2, std::move(N)) {}

BSplineCurve3D::BSplineCurve3D(vec K, field<vec>&& N)
    : BSpline(std::move(K), 3, std::move(N)) {}

BSplineCurve4D::BSplineCurve4D(vec K, field<vec>&& N)
    : BSpline(std::move(K), 4, std::move(N)) {}
