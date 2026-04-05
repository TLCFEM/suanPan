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

#include <vector>

template<std::floating_point T = double, unsigned BUCKET_SIZE = 1> class QuadTree {
    const BoundingBox<T> box;

    std::vector<Node2D<T>> nodes;

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
    }

    void quick_insert(Node2D<T>&& node) {
        const std::size_t a = node.x > box.center.x, b = node.y > box.center.y;
        children[2 * b + (a ^ b)].insert(std::move(node));
    }

public:
    explicit QuadTree(BoundingBox<T>&& in_box)
        : box(std::move(in_box)) { nodes.reserve(BUCKET_SIZE + 1); }

    void insert(Node2D<T>&& node) {
        if(!children.empty()) return quick_insert(std::move(node));

        nodes.emplace_back(std::move(node));
        if(nodes.size() <= BUCKET_SIZE) return;

        split();
        for(auto&& n : nodes) quick_insert(std::move(n));
        nodes.clear();
        nodes.shrink_to_fit();
    }
};

#endif
