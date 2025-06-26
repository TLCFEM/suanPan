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

#ifndef RIDDERS_HPP
#define RIDDERS_HPP

#include <suanPan.h>

/**
 * @brief Implements Ridders' method for finding the root of a function.
 *
 * Ridders' method is a root-finding algorithm that combines the bisection method
 * with an exponential interpolation to achieve quadratic convergence under ideal conditions.
 *
 * @tparam T A callable type that takes a double as input and returns a double.
 * @param func The function for which the root is to be found. It must be callable with a double argument.
 * @param x1 The first initial guess for the root.
 * @param f1 The value of the function at x1 (i.e., func(x1)).
 * @param x2 The second initial guess for the root.
 * @param f2 The value of the function at x2 (i.e., func(x2)).
 * @param tolerance The tolerance for the root-finding process.
 * @return The estimated root of the function within the specified tolerance.
 */
template<std::invocable<double> T> double ridders(const T& func, double x1, double f1, double x2, double f2, const double tolerance) {
    double target;

    auto counter = 0u;
    while(true) {
        counter += 2u;

        const auto x3 = .5 * (x1 + x2);
        const auto f3 = func(x3);
        if(std::fabs(f3) < tolerance || std::fabs(x2 - x1) < tolerance) {
            target = x3;
            break;
        }

        const auto dx = (x3 - x1) * f3 / std::sqrt(f3 * f3 - f1 * f2);

        const auto x4 = f1 > f2 ? x3 + dx : x3 - dx;
        const auto f4 = func(x4);
        if(std::fabs(f4) < tolerance) {
            target = x4;
            break;
        }

        // one end is x4
        // pick the other from x3, x2, x1
        if(std::signbit(f4) != std::signbit(f3)) {
            x1 = x3;
            f1 = f3;
        }
        else if(std::signbit(f4) != std::signbit(f2)) {
            x1 = x2;
            f1 = f2;
        }

        x2 = x4;
        f2 = f4;
    }

    suanpan_debug("Ridders' method initial guess {:.5E} with {} iterations.\n", target, counter);

    return target;
}

template<std::invocable<double> T> double ridders(const T& func, double x1, double x2, const double tolerance) { return ridders(func, x1, func(x1), x2, func(x2), tolerance); }

template<std::invocable<double> T> double ridders_guess(const T& func, double x1, double f1, double x2, double f2, const double tolerance) {
    while(std::signbit(f1) == std::signbit(f2)) {
        x1 = x2;
        f1 = f2;
        f2 = func(x2 *= 2.);
    }

    return ridders(func, x1, f1, x2, f2, tolerance);
}

template<std::invocable<double> T> double ridders_guess(const T& func, double x1, double x2, const double tolerance) { return ridders_guess(func, x1, func(x1), x2, func(x2), tolerance); }

#endif

//! @}
