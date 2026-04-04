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

#include <Toolbox/tree/pquadtree.hpp>
#include <iostream>
#include <random>
#include <ranges>

void test_mode() {
    std::mt19937 gen(42);
    std::uniform_real_distribution dis(-1000.0, 1000.0);

    std::vector<Node2D<>> points;
    points.reserve(30'000'000);
    for(auto i = 0; i < 30'000'000; ++i) points.push_back(Node2D<>{dis(gen), dis(gen)});

    auto ptr_view = points | std::views::transform([](const Node2D<>& p) { return &p; });

    QuadTree<double, 2> tree({{0.0, 0.0}, {1000.0, 1000.0}});
    tree.insert(ptr_view.begin(), ptr_view.end());

    std::cout << points.size() << std::endl;
}
