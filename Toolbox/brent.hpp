/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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

#ifndef BRENT_HPP
#define BRENT_HPP

#include <suanPan.h>

/**
 * @brief Implements Brent's method for finding a root of a function within a given interval.
 *
 * Brent's method combines the bisection method, the secant method, and inverse quadratic interpolation
 * to efficiently find a root of a function. It is robust and guarantees convergence as long as the
 * function changes sign over the interval [x1, x2].
 *
 * @tparam T The floating-point type (e.g., float, double, long double).
 * @tparam F The callable type representing the function to find the root of. Must accept a single
 *           argument of type T and return a value of type T.
 *
 * @param func The function for which the root is to be found. It must be continuous and change sign
 *             over the interval [x1, x2].
 * @param x1 The lower bound of the interval.
 * @param x2 The upper bound of the interval.
 * @param tol The desired tolerance for the root. The algorithm stops when the root is found to this
 *            precision.
 *
 * @return The approximate root of the function within the given interval and tolerance.
 *
 * @note The function assumes that the initial interval [x1, x2] contains a root (i.e., func(x1) and
 *       func(x2) have opposite signs). If this condition is not met, the behavior is undefined.
 */
template<std::floating_point T, std::invocable<T> F> T brent(F&& func, const T x1, const T x2, const T tol) {
    static constexpr auto eps = std::numeric_limits<T>::epsilon();
    auto a = x1, b = x2, c = x2;
    T fa = func(a), fb = func(b), fc = fb;

    T d{}, e{}, p{}, q{};

    auto counter = 0u;
    while(true) {
        ++counter;

        if(fb * fc > 0) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if(std::fabs(fc) < std::fabs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        const auto comp_tol = 2 * eps * std::fabs(b) + .5 * tol;
        const auto xm = .5 * (c - b);
        if(std::fabs(xm) <= comp_tol || std::fabs(fb) < tol) {
            suanpan_debug("Brent's method initial guess {:.5E} with {} iterations.\n", b, counter);
            return b;
        }
        if(std::fabs(e) >= comp_tol && std::fabs(fa) > std::fabs(fb)) {
            const auto s = fb / fa;
            if(a == c) {
                p = 2 * xm * s;
                q = 1 - s;
            }
            else {
                q = fa / fc;
                const auto r = fb / fc;
                p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1));
                q = (q - 1) * (r - 1) * (s - 1);
            }
            if(p > 0) q = -q;
            else p = std::fabs(p);
            if(2 * p < std::min(3 * xm * q - std::fabs(comp_tol * q), std::fabs(e * q))) {
                e = d;
                d = p / q;
            }
            else e = d = xm;
        }
        else e = d = xm;
        a = b;
        fa = fb;
        const auto dd = xm >= 0 ? comp_tol : -comp_tol;
        fb = func(b += std::fabs(d) > comp_tol ? d : dd);
    }
}

#endif

//! @}
