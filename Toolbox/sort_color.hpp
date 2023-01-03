/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include <set>
#include "container.h"
#include "utility.h"
#include "metis/metis.h"

template<typename T> auto sort_color_metis(suanpan::graph<T>& element_register, const int num_color, const char method) {
    wall_clock timer;
    timer.tic();

    const auto element_size = element_register.size();

    std::atomic num_edges = 0llu;
    suanpan::for_all(element_register, [&](const suanpan::set<T>& element) { num_edges += element.size(); });

    std::vector<idx_t> xadj;
    xadj.reserve(element_size + 1llu);
    xadj.emplace_back(0);
    std::vector<idx_t> adjncy;
    adjncy.reserve(num_edges);

    for(const auto& element : element_register) {
        adjncy.insert(adjncy.end(), element.begin(), element.end());
        xadj.emplace_back(xadj.back() + static_cast<idx_t>(element.size()));
    }

    auto nvtxs = static_cast<idx_t>(element_size);
    idx_t ncon = 1;
    auto nparts = static_cast<idx_t>(num_color);
    idx_t edgecut = 0;

    std::vector<int> part(nvtxs);

    idx_t* vsize = nullptr;
    real_t* tpwgts = nullptr;
    real_t* ubvec = nullptr;

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
#ifdef SUANPAN_DEBUG
    options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO | METIS_DBG_TIME;
#endif
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

    if('K' == method)
        METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(), nullptr, vsize, nullptr, &nparts, tpwgts, ubvec, options, &edgecut, part.data());
    else
        METIS_PartGraphRecursive(&nvtxs, &ncon, xadj.data(), adjncy.data(), nullptr, vsize, nullptr, &nparts, tpwgts, ubvec, options, &edgecut, part.data());

    SP_D("Coloring algorithm takes {:.5E} seconds.\n", timer.toc());

    return part;
}

template<typename T> std::vector<std::vector<T>> sort_color_wp(const suanpan::graph<T>& node_register) {
    wall_clock timer;
    timer.tic();

    const auto num_node = node_register.size();

    uvec weight(num_node, fill::none);

    suanpan_for(static_cast<size_t>(0), node_register.size(), [&](const size_t I) { weight(I) = node_register[I].size(); });

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

    SP_D("Coloring algorithm takes {:.5E} seconds.\n", timer.toc());

    return color_map;
}

template<typename T> std::vector<std::vector<T>> sort_color_mis(const suanpan::graph<T>& node_register) {
    wall_clock timer;
    timer.tic();

    const auto num_node = node_register.size();

    uvec weight(num_node, fill::none);

    suanpan_for(static_cast<size_t>(0), node_register.size(), [&](const size_t I) { weight(I) = node_register[I].size(); });

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
        }
        while(!degree_copy.empty());
    }

    for(auto& color : color_map) color.shrink_to_fit();
    color_map.shrink_to_fit();

    SP_D("Coloring algorithm takes {:.5E} seconds.\n", timer.toc());

    return color_map;
}

#endif

//! @}
