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

#ifndef QUADTREE_HPP
#define QUADTREE_HPP

#include "common.hpp"

#include <array>
#include <memory>
#include <variant>
#include <vector>

#ifdef SUANPAN_MT
#include <oneapi/tbb/parallel_for.h>
#endif

template<std::floating_point T, unsigned BUCKET_SIZE> class QuadTree;

template<std::floating_point T = double, unsigned BUCKET_SIZE = 1> struct RefNode2D : Node2D<T> {
    const QuadTree<T, BUCKET_SIZE>* tree{nullptr};
};

template<std::floating_point T = double, unsigned BUCKET_SIZE = 1> class QuadTree {
public:
    using node_t = RefNode2D<T, BUCKET_SIZE>;
    using node_ptr_t = node_t*;

private:
    using leaf_t = std::vector<node_ptr_t>;
    using subtree_t = std::array<std::unique_ptr<QuadTree>, 4>;

    const BoundingBox<T> box;

    std::variant<leaf_t, subtree_t> child;

    const QuadTree* parent{nullptr};

    void attach(const QuadTree* in_parent) { parent = in_parent; }

public:
    explicit QuadTree(BoundingBox<T>&& in_box)
        : box(std::move(in_box)) {}

    template<std::forward_iterator IT> void insert(IT begin, IT end) requires std::is_convertible_v<std::iter_value_t<IT>, node_ptr_t> { insert(leaf_t(begin, end)); }

    void insert(leaf_t&& in_child) {
        if(auto& nodes = std::get<leaf_t>(child = std::move(in_child)); nodes.size() <= BUCKET_SIZE) {
            for(auto&& node : nodes) node->tree = this;
            return;
        }

        std::array<leaf_t, 4> buckets;
        for(auto&& node : std::get<leaf_t>(child)) buckets[2 * (node->y > box.center.y) + (node->x > box.center.x)].push_back(node);
        for(auto i = 0; i < 4; ++i) buckets[i].shrink_to_fit();

        auto& subtrees = std::get<subtree_t>(child = subtree_t());

        const auto dx = T(.5) * box.dimension.x, dy = T(.5) * box.dimension.y;
        if(!buckets[0].empty()) {
            subtrees[0] = std::make_unique<QuadTree>(BoundingBox<T>{{box.center.x - dx, box.center.y - dy}, {dx, dy}});
            subtrees[0]->attach(this);
        }
        if(!buckets[1].empty()) {
            subtrees[1] = std::make_unique<QuadTree>(BoundingBox<T>{{box.center.x + dx, box.center.y - dy}, {dx, dy}});
            subtrees[1]->attach(this);
        }
        if(!buckets[2].empty()) {
            subtrees[2] = std::make_unique<QuadTree>(BoundingBox<T>{{box.center.x - dx, box.center.y + dy}, {dx, dy}});
            subtrees[2]->attach(this);
        }
        if(!buckets[3].empty()) {
            subtrees[3] = std::make_unique<QuadTree>(BoundingBox<T>{{box.center.x + dx, box.center.y + dy}, {dx, dy}});
            subtrees[3]->attach(this);
        }

#ifdef SUANPAN_MT
        tbb::parallel_for(0, 4, [&](const auto i) {
            if(subtrees[i]) subtrees[i]->insert(std::move(buckets[i]));
        });
#else
        for(auto i = 0; i < 4; ++i)
            if(subtrees[i]) subtrees[i]->insert(std::move(buckets[i]));
#endif
    }

    std::vector<const QuadTree*> overlap(const BoundingBox<T>& target) const {
        if(!box.overlap(target)) return {};

        if(std::holds_alternative<leaf_t>(child)) return {this};

        std::vector<const QuadTree*> result;
        for(auto&& subtree : std::get<subtree_t>(child)) {
            if(!subtree) continue;
            auto sub_result = subtree->overlap(target);
            result.insert(result.end(), sub_result.begin(), sub_result.end());
        }
        return result;
    }
};

#endif
