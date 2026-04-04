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

#include <oneapi/tbb/parallel_for.h>

template<std::floating_point T = double, unsigned BUCKET_SIZE = 1> class QuadTree {
    using node_ptr = const Node2D<T>*;
    using node_pool = std::vector<node_ptr>;

    const BoundingBox<T> box;

    node_pool nodes;

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

        std::array<node_pool, 4> buckets;
        for(auto&& node : nodes) {
            const std::size_t a = node->x > box.center.x, b = node->y > box.center.y;
            buckets[2 * b + (a ^ b)].push_back(node);
        }
        nodes.clear();
        nodes.shrink_to_fit();

        tbb::parallel_for(0, 4, [&](const auto i) { children[i].insert(std::move(buckets[i])); });
    }

    void insert(node_pool&& child) {
        nodes = std::move(child);

        if(nodes.size() > BUCKET_SIZE) split();
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
