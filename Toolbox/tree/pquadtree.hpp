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
    using leaf_t = std::vector<node_ptr>;
    using subtree_t = std::array<std::unique_ptr<PQuadTree>, 4>;

    const BoundingBox<T> box;

    std::variant<leaf_t, subtree_t> child;

    const PQuadTree* parent = nullptr;

    void attach(const PQuadTree* in_parent) { parent = in_parent; }

public:
    explicit PQuadTree(BoundingBox<T>&& in_box)
        : box(std::move(in_box)) {}

    template<std::forward_iterator IT> void insert(IT begin, IT end) requires std::is_convertible_v<std::iter_value_t<IT>, node_ptr> { insert(leaf_t(begin, end)); }

    void insert(leaf_t&& in_child) {
        if(std::get<leaf_t>(child = std::move(in_child)).size() <= BUCKET_SIZE) return;

        std::array<leaf_t, 4> buckets;
        for(auto&& node : std::get<leaf_t>(child)) buckets[2 * (node->y > box.center.y) + (node->x > box.center.x)].push_back(node);
        for(auto i = 0; i < 4; ++i) buckets[i].shrink_to_fit();

        auto& subtree = std::get<subtree_t>(child = subtree_t());

        const auto dx = T(.5) * box.dimension.x, dy = T(.5) * box.dimension.y;
        if(!buckets[0].empty()) {
            subtree[0] = std::make_unique<PQuadTree>(BoundingBox<T>{{box.center.x - dx, box.center.y - dy}, {dx, dy}});
            subtree[0]->attach(this);
        }
        if(!buckets[1].empty()) {
            subtree[1] = std::make_unique<PQuadTree>(BoundingBox<T>{{box.center.x + dx, box.center.y - dy}, {dx, dy}});
            subtree[1]->attach(this);
        }
        if(!buckets[2].empty()) {
            subtree[2] = std::make_unique<PQuadTree>(BoundingBox<T>{{box.center.x - dx, box.center.y + dy}, {dx, dy}});
            subtree[2]->attach(this);
        }
        if(!buckets[3].empty()) {
            subtree[3] = std::make_unique<PQuadTree>(BoundingBox<T>{{box.center.x + dx, box.center.y + dy}, {dx, dy}});
            subtree[3]->attach(this);
        }

#ifdef SUANPAN_MT
        tbb::parallel_for(0, 4, [&](const auto i) {
            if(subtree[i]) subtree[i]->insert(std::move(buckets[i]));
        });
#else
        for(auto i = 0; i < 4; ++i)
            if(subtree[i]) subtree[i]->insert(std::move(buckets[i]));
#endif
    }
};

#endif
