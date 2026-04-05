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

#include <variant>

#ifdef SUANPAN_MT
#include <oneapi/tbb/parallel_for.h>
#endif

template<std::floating_point T = double, unsigned BUCKET_SIZE = 1> class PQuadTree {
    using node_ptr = const Node2D<T>*;
    using node_pool = std::vector<node_ptr>;

    const BoundingBox<T> box;

    std::variant<node_pool, std::vector<PQuadTree>> child;

    const PQuadTree* parent = nullptr;

    void attach(const PQuadTree* in_parent) { parent = in_parent; }

public:
    explicit PQuadTree(BoundingBox<T>&& in_box)
        : box(std::move(in_box)) {}

    template<std::forward_iterator IT> void insert(IT begin, IT end) requires std::is_convertible_v<std::iter_value_t<IT>, node_ptr> { insert(node_pool(begin, end)); }

    void insert(node_pool&& in_child) {
        if(std::get<node_pool>(child = std::move(in_child)).size() <= BUCKET_SIZE) return;

        std::array<node_pool, 4> buckets;
        for(auto&& node : std::get<node_pool>(child)) buckets[2 * (node->y > box.center.y) + (node->x > box.center.x)].push_back(node);
        for(auto i = 0; i < 4; ++i) buckets[i].shrink_to_fit();

        auto& subtree = std::get<std::vector<PQuadTree>>(child = std::vector<PQuadTree>());
        subtree.reserve(4);

        const auto dx = T(.5) * box.dimension.x, dy = T(.5) * box.dimension.y;
        subtree.emplace_back(BoundingBox<T>{{box.center.x - dx, box.center.y - dy}, {dx, dy}}).attach(this);
        subtree.emplace_back(BoundingBox<T>{{box.center.x + dx, box.center.y - dy}, {dx, dy}}).attach(this);
        subtree.emplace_back(BoundingBox<T>{{box.center.x - dx, box.center.y + dy}, {dx, dy}}).attach(this);
        subtree.emplace_back(BoundingBox<T>{{box.center.x + dx, box.center.y + dy}, {dx, dy}}).attach(this);

#ifdef SUANPAN_MT
        tbb::parallel_for(0, 4, [&](const auto i) { subtree[i].insert(std::move(buckets[i])); });
#else
        for(auto i = 0; i < 4; ++i) subtree[i].insert(std::move(buckets[i]));
#endif
    }
};

#endif
