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

uvec sort_rcm(const std::vector<uvec>& A, const uvec& E) {
#ifdef SUANPAN_DEBUG
    wall_clock TM;
    TM.tic();
#endif

    const auto S = E.n_elem;

    uvec G = sort_index(E);
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

#ifdef SUANPAN_DEBUG
    suanpan_debug("RCM algorithm takes %.5E seconds.\n", TM.toc());
#endif

    return R;
}
