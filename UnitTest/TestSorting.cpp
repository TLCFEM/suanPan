/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include <Toolbox/sort_rcm.h>
#include "CatchHeader.h"

TEST_CASE("RCM", "[Utility.Sorting]") {
    constexpr auto N = 200;

    for(auto I = 0; I < 4; ++I) {
        sp_mat B = sprandu(N, N, .01);
        sp_mat A = B + B.t();

        auto C = mat(A);

        BENCHMARK("RCM Dense") { return sort_rcm(C); };

        triplet_form<double, uword> triplet_mat(N, N);

        for(auto J = A.begin(); J != A.end(); ++J) triplet_mat.at(J.row(), J.col()) = *J;

        BENCHMARK("RCM Triplet") { return sort_rcm(triplet_mat); };

        csc_form<double, uword> csc_mat(triplet_mat);

        BENCHMARK("RCM CSC") { return sort_rcm(csc_mat); };
    }
}
