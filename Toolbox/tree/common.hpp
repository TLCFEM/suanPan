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

#ifndef TREE_COMMON_HPP
#define TREE_COMMON_HPP

#include <cmath>
#include <concepts>
#include <cstdint>

template<std::floating_point T = double> struct Vector2D {
    const T x, y;
};

template<std::floating_point T = double> struct Node2D : Vector2D<T> {
    const std::uint64_t id{0};
};

template<std::floating_point T = double> struct BoundingBox {
    const Vector2D<T> center, dimension;

    bool overlap(const BoundingBox& other) const { return std::fabs(center.x - other.center.x) < (dimension.x + other.dimension.x) && std::fabs(center.y - other.center.y) < (dimension.y + other.dimension.y); }
};

#endif
