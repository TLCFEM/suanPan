/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#ifndef SUBLOADINGUTIL_H
#define SUBLOADINGUTIL_H

#include <array>
#include <cmath>
#include <utility>

struct SubloadingBound {
    const double initial, linear, saturation, rate;

    [[nodiscard]] std::pair<double, double> operator()(const double q, const bool check_positive) const {
        const auto exp_term = saturation * std::exp(-rate * q);
        const auto y = initial + saturation + linear * q - exp_term;
        return y < 0. && check_positive ? std::make_pair(0., 0.) : std::make_pair(y, linear + rate * exp_term);
    }
};

class SubloadingSaturation {
    const double rate, bound;

public:
    SubloadingSaturation(const double R, const double B)
        : rate(R)
        , bound(B) {}

    [[nodiscard]] double r() const { return rate; }
    [[nodiscard]] double b() const { return bound; }
    [[nodiscard]] double rb() const { return r() * b(); }
};

class SubloadingBase {
    static constexpr double z_bound = 1E-15;
    inline static const double rate_bound = -std::log(z_bound);

protected:
    static std::array<double, 2> yield_ratio(const double z) {
        if(z < z_bound) return {rate_bound, 0.};

        return {-log(z), -1. / z};
    }
};

#endif
