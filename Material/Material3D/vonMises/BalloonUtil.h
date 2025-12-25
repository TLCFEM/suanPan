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

#ifndef BALLOONUTIL_H
#define BALLOONUTIL_H

#include <Material/Material.h>
#include <numeric> // std::accumulate

class BalloonBuffer {
    std::vector<double> buffer;
    std::size_t head;

    [[nodiscard]] auto max() const { return *std::ranges::max_element(buffer); }
    [[nodiscard]] auto min() const { return *std::ranges::min_element(buffer); }
    [[nodiscard]] auto mean() const { return std::accumulate(buffer.cbegin(), buffer.cend(), 0.) / static_cast<double>(buffer.size()); }

public:
    enum class Type : std::uint8_t {
        MINIMUM,
        MAXIMUM,
        MEAN
    };

    explicit BalloonBuffer(const std::size_t size)
        : buffer(std::max(std::size_t{1}, size), 0.)
        , head(0) {}

    auto& enqueue(const double value) {
        buffer[head] = value;
        head = (head + 1) % buffer.size();
        return *this;
    }

    auto zeros() {
        std::ranges::fill(buffer, 0.);
        head = 0;
    }

    [[nodiscard]] auto operator()(const Type memory_type) const {
        switch(memory_type) {
        case Type::MINIMUM:
            return min();
        case Type::MAXIMUM:
            return max();
        case Type::MEAN:
        default:
            return mean();
        }
    }
};

class BalloonBound {
    const double initial, linear, saturation, rate;

public:
    BalloonBound(const double I, const double K, const double S, const double M)
        : initial(I)
        , linear(K)
        , saturation(S)
        , rate(M) {}

    [[nodiscard]] std::pair<double, double> operator()(const double q, const bool check_positive) const {
        const auto exp_term = saturation * std::exp(-rate * q);
        const auto y = initial + saturation + linear * q - exp_term;
        return y < 0. && check_positive ? std::make_pair(0., 0.) : std::make_pair(y, linear + rate * exp_term);
    }
};

class BalloonSaturation {
    const double rate, bound;

public:
    BalloonSaturation(const double R, const double B)
        : rate(R)
        , bound(B) {}

    [[nodiscard]] double a() const { return (rate > 0. ? b() : 1.) * bound; }
    [[nodiscard]] double b() const { return rate; }
};

class BalloonBase {
    static constexpr double z_bound = 1E-15;
    inline static const double rate_bound = -std::log(z_bound);

protected:
    static pod2 yield_ratio(const double z) {
        if(z < z_bound) return {rate_bound, 0.};

        return {-log(z), -1. / z};
    }
};

#endif
