/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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
/**
 * @fn sort_color
 * @brief A color sorting algorithm.
 * @author tlc
 * @date 22/02/2022
 * @version 0.1.2
 * @file sort_color.hpp
 * @addtogroup Utility
 * @{
 */

// ReSharper disable IdentifierTypo
#ifndef SORT_COLOR_HPP
#define SORT_COLOR_HPP

#include "container.h"
#include "utility.h"

#include <set>

template<typename T> std::vector<std::vector<T>> sort_color_wp(const suanpan::graph<T>& node_register) {
    wall_clock timer;
    timer.tic();

    const auto num_node = node_register.size();

    uvec weight(num_node, fill::none);

    suanpan::for_each(node_register.size(), [&](const size_t I) { weight(I) = node_register[I].size(); });

    auto comparator = [&](const T A, const T B) { return weight[A] > weight[B]; };

    std::multiset<T, decltype(comparator)> degree(comparator);
    for(T I = 0; I < num_node; ++I) degree.insert(I);

    std::vector<std::vector<T>> color_map;

    while(!degree.empty()) {
        color_map.emplace_back();
        auto& color = color_map.back();

        color.emplace_back(*degree.begin());
        for(auto node = degree.erase(degree.begin()); node != degree.end();) {
            auto flag = false;
            for(const auto colored : color) {
                if(colored == *node) continue;
                if(node_register[colored].contains(*node)) {
                    flag = true;
                    break;
                }
            }
            if(flag) ++node;
            else {
                color.emplace_back(*node);
                node = degree.erase(node);
            }
        }
    }

    for(auto& color : color_map) color.shrink_to_fit();
    color_map.shrink_to_fit();

    suanpan_debug("Coloring algorithm takes {:.5E} seconds.\n", timer.toc());

    return color_map;
}

template<typename T> std::vector<std::vector<T>> sort_color_mis(const suanpan::graph<T>& node_register) {
    wall_clock timer;
    timer.tic();

    const auto num_node = node_register.size();

    uvec weight(num_node, fill::none);

    suanpan::for_each(node_register.size(), [&](const size_t I) { weight(I) = node_register[I].size(); });

    uword counter = num_node;
    for(const auto I : sort_index(weight, "descend").eval()) weight[I] = --counter;

    auto comparator = [&](const T A, const T B) { return weight[A] > weight[B]; };

    std::multiset<T, decltype(comparator)> degree(comparator);
    for(T I = 0; I < num_node; ++I) degree.insert(I);

    std::vector<std::vector<T>> color_map;

    while(!degree.empty()) {
        color_map.emplace_back();
        auto& color = color_map.back();

        auto degree_copy = degree;
        do {
            const auto node_to_color = *degree_copy.begin();
            color.emplace_back(node_to_color);
            degree.erase(node_to_color);
            for(auto neighbour : node_register[node_to_color]) degree_copy.erase(neighbour);
        } while(!degree_copy.empty());
    }

    for(auto& color : color_map) color.shrink_to_fit();
    color_map.shrink_to_fit();

    suanpan_debug("Coloring algorithm takes {:.5E} seconds.\n", timer.toc());

    return color_map;
}

#endif

//! @}
