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

#include "sort_rcm.h"
#include <unordered_set>

/**
 * \brief Find a pair of pseudo-peripheral vertices
 *
 * See the Gibbs-Poole-Stockmeyer (GPS) algorithm.
 *
 * https://doi.org/10.1007/978-3-031-25820-6
 *
 * \param A adjacency list
 * \param E vertex degree list
 * \return a list of vertices beginning with the pseudo-peripheral vertex
 */
uvec peripheral(const std::vector<uvec>& A, const uvec& E) {
    uvec order = sort_index(E);

    std::size_t depth = 0;
    auto root = order(0);

    while(true) {
        std::vector mask(E.n_elem, false);

        // ReSharper disable once CppTemplateArgumentsCanBeDeduced
        std::unordered_set<decltype(root)> parent_set{root};
        mask[root] = true;
        std::size_t current_depth = 1;

        const auto populate = [&] {
            decltype(parent_set) child_set;

            for(const auto parent : parent_set)
                for(const auto child : A[parent])
                    if(!mask[child]) {
                        child_set.insert(child);
                        mask[child] = true;
                    }

            return child_set;
        };

        // loop over all vertices and construct the level set
        while(true) {
            auto child_set = populate();
            if(child_set.empty()) break;
            parent_set = std::move(child_set);
            ++current_depth;
        }

        if(current_depth <= depth) break;

        root = *parent_set.begin();
        depth = current_depth;
    }

    std::swap(*order.begin(), *std::find(order.begin(), order.end(), root));

    return order;
}

uvec sort_rcm(const std::vector<uvec>& A, const uvec& E) {
    wall_clock TM;
    TM.tic();

    const auto S = E.n_elem;

    const auto G = peripheral(A, E);
    uvec R(S, fill::zeros);
    std::vector M(S, false);

    uword IDXA = 0, IDXB = S - 1, IDXC = S - 1;

    while(true) {
        if(IDXB == IDXC) {
            while(IDXA < S && M[G(IDXA)]) ++IDXA;
            if(IDXA == S) break;
            M[R(IDXC--) = G(IDXA++)] = true;
        }
        for(const auto& IDX : A[R(IDXB--)]) if(!M[IDX]) M[R(IDXC--) = IDX] = true;
    }

    suanpan_debug("RCM algorithm takes {:.5E} seconds.\n", TM.toc());

    return R;
}

uvec sort_rcm(const std::vector<suanpan::unordered_set<uword>>& adjacency) {
    const auto dof_size = adjacency.size();

    // count number of degree
    uvec num_degree(dof_size, fill::none);
    suanpan::for_each(static_cast<size_t>(0), dof_size, [&](const size_t I) { num_degree(I) = adjacency[I].size(); });

    // sort each column according to its degree
    std::vector<uvec> adjacency_sorted(dof_size);
    suanpan::for_each(static_cast<size_t>(0), dof_size, [&](const size_t I) {
        adjacency_sorted[I] = to_uvec(adjacency[I]);
        auto& t_vec = adjacency_sorted[I];
        std::sort(t_vec.begin(), t_vec.end(), [&](const uword a, const uword b) { return num_degree(a) < num_degree(b); });
    });

    return sort_rcm(adjacency_sorted, num_degree);
}
