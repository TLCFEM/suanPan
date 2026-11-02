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

#include "CatchHeader.h"

#include <Toolbox/sort_rcm.h>

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

TEST_CASE("RCM Is A Valid Permutation", "[Utility.Sorting]") {
    constexpr uword N = 100;

    for(auto I = 0; I < 10; ++I) {
        sp_mat B = sprandu(N, N, .02);
        sp_mat A = B + B.t();

        auto R = sort_rcm(A);

        REQUIRE(N == R.n_elem);

        // R must be a bijection of {0,...,N-1}: every value appears exactly once
        uvec sorted = sort(R);
        for(uword J = 0; J < N; ++J) REQUIRE(J == sorted(J));
    }
}

TEST_CASE("RCM Of Dense Graph Is A Valid Permutation", "[Utility.Sorting]") {
    // a fully connected small graph - only checks the permutation validity property
    constexpr uword N = 30;
    sp_mat A(N, N);
    for(uword I = 0; I < N; ++I)
        for(uword J = 0; J < N; ++J)
            if(I != J) A(I, J) = 1.;

    const auto R = sort_rcm(A);
    REQUIRE(N == R.n_elem);
    uvec sorted = sort(R);
    for(uword J = 0; J < N; ++J) REQUIRE(J == sorted(J));
}
