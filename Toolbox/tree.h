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

#ifndef TREE_H
#define TREE_H

#include <cmath>
#include <vector>

struct Vector2D {
    double x, y;
};

struct Node2D : public Vector2D {
    unsigned id{0};
};

struct BoundingBox {
    Vector2D center, dimension;
};

template<unsigned BUCKET_SIZE = 1> class QuadTree {
    BoundingBox box;

    std::vector<Node2D> nodes;
    std::vector<QuadTree> children;

    const QuadTree* parent = nullptr;

    void attach(const QuadTree* in_parent) { parent = in_parent; }

    bool cover(const Node2D& node) const { return std::fabs(node.x - box.center.x) <= box.dimension.x && std::fabs(node.y - box.center.y) <= box.dimension.y; }

    void split() {
        const auto dx = .5 * box.dimension.x, dy = .5 * box.dimension.y;
        children.clear();
        children.reserve(4);
        children.emplace_back(BoundingBox{{box.center.x - dx, box.center.y - dy}, {dx, dy}}).attach(this);
        children.emplace_back(BoundingBox{{box.center.x + dx, box.center.y - dy}, {dx, dy}}).attach(this);
        children.emplace_back(BoundingBox{{box.center.x + dx, box.center.y + dy}, {dx, dy}}).attach(this);
        children.emplace_back(BoundingBox{{box.center.x - dx, box.center.y + dy}, {dx, dy}}).attach(this);
    }

public:
    QuadTree(BoundingBox&& in_box)
        : box(std::move(in_box)) {}

    void insert(Node2D&& node) {
        if(!cover(node)) return;

        if(children.empty()) {
            nodes.emplace_back(std::move(node));
            if(nodes.size() > BUCKET_SIZE) {
                split();
                for(auto&& n : nodes) insert(std::move(n));
                nodes.clear();
            }
        }
        else if(node.y <= box.center.y) children[node.x <= box.center.x ? 0 : 1].insert(std::move(node));
        else children[node.x <= box.center.x ? 3 : 2].insert(std::move(node));
    }
};

#endif
