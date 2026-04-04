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

#ifndef PQUADTREE_HPP
#define PQUADTREE_HPP

#include "common.hpp"

#include <oneapi/tbb/concurrent_vector.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_for_each.h>

template<std::floating_point T = double, unsigned BUCKET_SIZE = 1> class QuadTree {
    using node_ptr = const Node2D<T>*;

    const BoundingBox<T> box;

    tbb::concurrent_vector<node_ptr> nodes;

    std::vector<QuadTree> children;

    const QuadTree* parent = nullptr;

    void attach(const QuadTree* in_parent) { parent = in_parent; }

    void split() {
        const auto dx = T(.5) * box.dimension.x, dy = T(.5) * box.dimension.y;
        children.clear();
        children.reserve(4);
        children.emplace_back(BoundingBox<T>{{box.center.x - dx, box.center.y - dy}, {dx, dy}}).attach(this);
        children.emplace_back(BoundingBox<T>{{box.center.x + dx, box.center.y - dy}, {dx, dy}}).attach(this);
        children.emplace_back(BoundingBox<T>{{box.center.x + dx, box.center.y + dy}, {dx, dy}}).attach(this);
        children.emplace_back(BoundingBox<T>{{box.center.x - dx, box.center.y + dy}, {dx, dy}}).attach(this);

        std::array<decltype(nodes), 4> buckets;
        tbb::parallel_for_each(nodes.cbegin(), nodes.cend(), [&](const auto& node) {
            const std::size_t a = node->x > box.center.x, b = node->y > box.center.y;
            buckets[2 * b + (a ^ b)].push_back(node);
        });
        nodes.clear();

        tbb::parallel_for(0, 4, [&](const auto i) {
            children[i].insert(buckets[i].cbegin(), buckets[i].cend());
            buckets[i].clear();
        });
    }

public:
    explicit QuadTree(BoundingBox<T>&& in_box)
        : box(std::move(in_box)) {}

    template<std::forward_iterator IT> void insert(IT begin, IT end) requires std::is_convertible_v<std::iter_value_t<IT>, node_ptr> {
        nodes.assign(begin, end);

        if(nodes.size() > BUCKET_SIZE) split();
    }
};

#endif
