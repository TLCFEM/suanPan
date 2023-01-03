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
 * @fn sort_rcm
 * @brief A renumber function using RCM algorithm.
 *
 * The function takes both mat and sp_mat.
 *
 * Example Usage:
 *
 * ```cpp
 *     sp_mat test_rcm=sprandn(100000, 100000, 0.00005);
 *     auto R = RCM(test_rcm + test_rcm.t());
 * ```
 *
 * R gives the new numbering order of the original symmetric matrix.
 *
 * @author tlc
 * @date 02/08/2017
 * @version 0.1.2
 * @file sort_rcm.h
 * @addtogroup Utility
 * @{
 */

#ifndef RCM_H
#define RCM_H

#include <Domain/MetaMat/triplet_form.hpp>
#include <Domain/MetaMat/csc_form.hpp>

uvec sort_rcm(const std::vector<uvec>&, const uvec&);

template<typename eT> uvec sort_rcm(const SpMat<eT>& MEAT) {
    suanpan_debug([&] {if(!MEAT.is_square()) throw logic_error("RCM() can only be applied to square matrix.\n"); });

    wall_clock TM;
    TM.tic();

    //! Get the size of the square matrix.
    auto S = MEAT.n_cols;

    //! Collect the number of degree of each node.
    uvec E(S, fill::none);
    suanpan_for(0llu, S, [&](const uword I) { E(I) = MEAT.col(I).n_nonzero; });

    std::vector<uvec> A(S);
    suanpan_for(0llu, S, [&](const uword K) {
        unsigned J = 0;
        uvec IDX(E(K));
        for(auto L = MEAT.begin_col(K); L != MEAT.end_col(K); ++L) IDX(J++) = L.row();
        A[K] = IDX(sort_index(E(IDX)));
    });

    //! Get the indices array in increasing order of degree.
    uvec G = sort_index(E);

    //! Now define the mask vector to indicate if one node is numbered or not.
    std::vector<bool> M(S, false);
    //! Define the new order vector.
    uvec R(S, fill::zeros);

    //! Preparation is ready.
    //! The G stores all vertices increasing order of degree.
    //! The adjacency stores the neighbor vertices for each vertex sorted
    //! according to the number of degree. Now start to loop over all nodes.

    //! Define one position indicator IDXA for numbered vertices and two other
    //! position indicator IDXB and IDXC for looping.
    uword IDXA = 0, IDXB = S - 1, IDXC = S - 1;

    //! While the sum of mask does not equal to the size of matrix, there are
    //! vertices not numbered in the system, the algorithm should continue.

    //! IDXC will always point to the end of the vector. For any time when
    //! IDXB==IDXC, both indicators reach the end of reorder vector, i.e., there
    //! is no sub level any more (forms an independent subset). The graph cannot
    //! grow any more, then get the vertex with minimum degree as the new start
    //! point.
    while(IDXA < S) {
        if(IDXB == IDXC) {
            //! IDXA should be less that S so that G(IDXA) does not overflow.
            while(IDXA < S && M[G(IDXA)]) ++IDXA;
            //! Once IDXA hits S, there is no unnumbered vertex in the graph. Quit the
            //! loop.
            if(IDXA == S) break;
            //! Push in first unnumbered element in the list.
            R(IDXC--) = G(IDXA);
            //! Label it as renumbered and move IDXA to next position.
            M[G(IDXA++)] = true;
        }
        //! Now we at least has one root, which is indicated by the indicator IDXB,
        //! in our graph, push in all children into the vector. As they are already
        //! sorted, we can simply push in. When the loop is finished, move IDXB to
        //! next position, which may have another root or the children of current
        //! root.
        for(const auto& IDX : A[R(IDXB--)])
            if(!M[IDX]) M[R(IDXC--) = IDX] = true;
    }

    SP_D("RCM algorithm takes {:.5E} seconds.\n", TM.toc());

    return R;
}

template<typename eT> uvec sort_rcm(const Mat<eT>& MEAT) { return sort_rcm(SpMat<eT>(MEAT)); }

template<typename dt> uvec sort_rcm(const csc_form<dt, uword>& csc_mat) {
    //! Get the size of the square matrix.
    auto S = csc_mat.n_cols;

    //! Collect the number of degree of each node.
    uvec E(S, fill::none);
    suanpan_for(0llu, S, [&](const uword I) { E(I) = csc_mat.col(I + 1) - csc_mat.col(I); });
    std::vector<uvec> A(S);
    suanpan_for(0llu, S, [&](const uword I) {
        const uvec IDX(csc_mat.row_mem() + csc_mat.col(I), E(I));
        A[I] = IDX(sort_index(E(IDX)));
    });

    return sort_rcm(A, E);
}

template<typename dt, typename it> uvec sort_rcm(triplet_form<dt, it>& triplet_mat) {
    csc_form<dt, uword> csc_mat(triplet_mat);

    return sort_rcm(csc_mat);
}

#endif

//! @}
