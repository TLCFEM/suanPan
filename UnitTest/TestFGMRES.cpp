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
/**
 * @file TestFGMRES.cpp
 * @brief Unit tests for the oneMKL FGMRES wrapper.
 *
 * @author tlc
 * @date 24/11/2026
 * @version 0.1.0
 * @addtogroup Utility
 * @{
 */

#include "CatchHeader.h"

#include <Toolbox/fgmres.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

    struct DenseKernel {
        const mat* matrix;
        uword n_rows;

        [[nodiscard]] mat operator*(const mat& x) const { return *matrix * x; }
    };

    // Always returns a matrix with the wrong number of columns.
    struct MalformedKernel {
        uword n_rows;

        [[nodiscard]] mat operator*(const mat&) const { return mat(n_rows, 2, arma::fill::ones); }
    };

    mat make_diagonally_dominant(const uword n) {
        mat a(n, n, arma::fill::zeros);
        for(uword i = 0; i < n; ++i) {
            a(i, i) = 10. + static_cast<double>(i);
            if(i + 1 < n) {
                a(i, i + 1) = -1.;
                a(i + 1, i) = -1.;
            }
            if(i + 2 < n) {
                a(i, i + 2) = .5;
                a(i + 2, i) = .5;
            }
        }
        return a;
    }

    mat make_nonsymmetric(const uword n) {
        mat a(n, n, arma::fill::zeros);
        for(uword i = 0; i < n; ++i) {
            a(i, i) = 4. + static_cast<double>(i);
            if(i + 1 < n) a(i, i + 1) = -2.;
            if(i > 0) a(i, i - 1) = -.5;
        }
        return a;
    }

    vec make_reference(const uword n) {
        vec out(n);
        for(uword i = 0; i < n; ++i) out[i] = static_cast<double>(i + 1);
        return out;
    }

    double relative_residual(const mat& a, const vec& x, const vec& b) { return arma::norm(a * x - b, 2) / arma::norm(b, 2); }

    double max_difference(const vec& lhs, const vec& rhs) {
        REQUIRE(lhs.n_elem == rhs.n_elem);
        double difference = 0.;
        for(uword i = 0; i < lhs.n_elem; ++i) difference = std::max(difference, std::abs(lhs[i] - rhs[i]));
        return difference;
    }

} // namespace

TEST_CASE("FGMRES convergence", "[Utility.FGMRES]") {
    constexpr uword n = 8;
    const auto a = make_diagonally_dominant(n);
    const auto expected = make_reference(n);
    vec b = a * expected;
    const DenseKernel kernel{&a, n};
    const vec diagonal = a.diag();

    SECTION("converges on a diagonally dominant system") {
        vec x(n, arma::fill::zeros);
        MKL_INT iterations = -1;

        REQUIRE(0 == fgmres_solve(kernel, diagonal, x.memptr(), b.memptr(), 1E-12, &iterations));
        REQUIRE(iterations > 0);
        REQUIRE(iterations <= static_cast<MKL_INT>(n));
        REQUIRE_THAT(relative_residual(a, x, b), Catch::Matchers::WithinAbs(0., 1E-10));
        REQUIRE_THAT(max_difference(x, expected), Catch::Matchers::WithinAbs(0., 1E-8));
    }

    SECTION("converges on a non-symmetric system") {
        const auto asymmetric = make_nonsymmetric(n);
        const auto reference = make_reference(n);
        vec rhs = asymmetric * reference;
        const DenseKernel asymmetric_kernel{&asymmetric, n};
        const vec asymmetric_diagonal = asymmetric.diag();

        vec x(n, arma::fill::zeros);
        MKL_INT iterations = -1;

        REQUIRE(0 == fgmres_solve(asymmetric_kernel, asymmetric_diagonal, x.memptr(), rhs.memptr(), 1E-12, &iterations));
        REQUIRE(iterations > 0);
        REQUIRE_THAT(relative_residual(asymmetric, x, rhs), Catch::Matchers::WithinAbs(0., 1E-10));
        REQUIRE_THAT(max_difference(x, reference), Catch::Matchers::WithinAbs(0., 1E-8));
    }

    SECTION("honours a nonzero initial guess") {
        const auto initial = .9 * expected;
        vec x = initial;
        MKL_INT iterations = -1;

        REQUIRE(0 == fgmres_solve(kernel, diagonal, x.memptr(), b.memptr(), 1E-12, &iterations));
        REQUIRE(iterations > 0);
        REQUIRE_THAT(relative_residual(a, x, b), Catch::Matchers::WithinAbs(0., 1E-10));
    }

    SECTION("does not modify the right-hand side") {
        const vec original_b = b;
        vec x(n, arma::fill::zeros);

        REQUIRE(0 == fgmres_solve(kernel, diagonal, x.memptr(), b.memptr(), 1E-12));
        REQUIRE_THAT(max_difference(b, original_b), Catch::Matchers::WithinAbs(0., 0.));
    }

    SECTION("iteration count is optional") {
        vec x(n, arma::fill::zeros);

        REQUIRE(0 == fgmres_solve(kernel, diagonal, x.memptr(), b.memptr(), 1E-12));
        REQUIRE_THAT(relative_residual(a, x, b), Catch::Matchers::WithinAbs(0., 1E-10));
    }

    SECTION("handles several right-hand sides independently") {
        for(uword i = 0; i < 3; ++i) {
            const vec scaled = static_cast<double>(i + 2) * expected;
            vec rhs = a * scaled;
            vec x(n, arma::fill::zeros);

            REQUIRE(0 == fgmres_solve(kernel, diagonal, x.memptr(), rhs.memptr(), 1E-12));
            REQUIRE_THAT(max_difference(x, scaled), Catch::Matchers::WithinAbs(0., 1E-8));
        }
    }

    SECTION("handles a system larger than the default Krylov dimension") {
        constexpr uword large_n = 160;
        const auto large_a = make_diagonally_dominant(large_n);
        const auto large_expected = make_reference(large_n);
        vec large_b = large_a * large_expected;
        const DenseKernel large_kernel{&large_a, large_n};
        const vec large_diagonal = large_a.diag();

        vec x(large_n, arma::fill::zeros);
        MKL_INT iterations = -1;

        REQUIRE(0 == fgmres_solve(large_kernel, large_diagonal, x.memptr(), large_b.memptr(), 1E-12, &iterations));
        REQUIRE(iterations > 0);
        REQUIRE_THAT(relative_residual(large_a, x, large_b), Catch::Matchers::WithinAbs(0., 1E-10));
        REQUIRE_THAT(max_difference(x, large_expected), Catch::Matchers::WithinAbs(0., 1E-8));
    }
}

TEST_CASE("FGMRES input validation", "[Utility.FGMRES]") {
    constexpr uword n = 8;
    const auto a = make_diagonally_dominant(n);
    vec b = a * make_reference(n);
    const DenseKernel kernel{&a, n};
    const vec diagonal = a.diag();
    vec x(n, arma::fill::zeros);

    SECTION("rejects a zero-sized problem") {
        const mat empty(0, 0);
        const DenseKernel empty_kernel{&empty, 0};
        double dummy = 0.;

        REQUIRE(-1 == fgmres_solve(empty_kernel, vec(), &dummy, &dummy, 1E-12));
    }

    SECTION("rejects null vectors") {
        REQUIRE(-1 == fgmres_solve(kernel, diagonal, nullptr, b.memptr(), 1E-12));
        REQUIRE(-1 == fgmres_solve(kernel, diagonal, x.memptr(), nullptr, 1E-12));
    }

    SECTION("rejects invalid tolerances") {
        REQUIRE(-1 == fgmres_solve(kernel, diagonal, x.memptr(), b.memptr(), 0.));
        REQUIRE(-1 == fgmres_solve(kernel, diagonal, x.memptr(), b.memptr(), -1E-12));
        REQUIRE(-1 == fgmres_solve(kernel, diagonal, x.memptr(), b.memptr(), std::numeric_limits<double>::quiet_NaN()));
        REQUIRE(-1 == fgmres_solve(kernel, diagonal, x.memptr(), b.memptr(), std::numeric_limits<double>::infinity()));
    }

    SECTION("rejects a mis-sized preconditioner") {
        REQUIRE(-1 == fgmres_solve(kernel, vec(n + 1), x.memptr(), b.memptr(), 1E-12));
        REQUIRE(-1 == fgmres_solve(kernel, vec(n - 1), x.memptr(), b.memptr(), 1E-12));
        REQUIRE(-1 == fgmres_solve(kernel, vec(), x.memptr(), b.memptr(), 1E-12));
    }

    SECTION("rejects zero and non-finite preconditioner entries") {
        vec zero_diagonal = diagonal;
        zero_diagonal[0] = 0.;
        REQUIRE(-1 == fgmres_solve(kernel, zero_diagonal, x.memptr(), b.memptr(), 1E-12));

        vec nan_diagonal = diagonal;
        nan_diagonal[1] = std::numeric_limits<double>::quiet_NaN();
        REQUIRE(-1 == fgmres_solve(kernel, nan_diagonal, x.memptr(), b.memptr(), 1E-12));
    }
}

TEST_CASE("FGMRES failure reporting", "[Utility.FGMRES]") {
    constexpr uword n = 8;
    const auto a = make_diagonally_dominant(n);
    const auto expected = make_reference(n);
    vec b = a * expected;
    const DenseKernel kernel{&a, n};
    const vec diagonal = a.diag();

    SECTION("reports non-convergence for an unattainable tolerance") {
        vec x(n, arma::fill::zeros);
        MKL_INT iterations = -1;

        REQUIRE(-7 == fgmres_solve(kernel, diagonal, x.memptr(), b.memptr(), 1E-30, &iterations));
        REQUIRE(iterations > 0);
        // The best approximation is still returned to the caller.
        REQUIRE_THAT(relative_residual(a, x, b), Catch::Matchers::WithinAbs(0., 1E-8));
    }

    SECTION("reports a malformed matrix-vector product") {
        const MalformedKernel malformed{n};
        vec x(n, arma::fill::zeros);

        REQUIRE(-4 == fgmres_solve(malformed, diagonal, x.memptr(), b.memptr(), 1E-12));
    }

    SECTION("handles a singular system with a consistent right-hand side") {
        // Two identical rows, but a nonzero diagonal so the preconditioner is valid.
        mat singular(n, n, arma::fill::eye);
        singular(0, 1) = 1.;
        singular.row(1) = singular.row(0);
        vec rhs(n, arma::fill::zeros);
        const DenseKernel singular_kernel{&singular, n};
        const vec singular_diagonal = singular.diag();

        vec x(n, arma::fill::ones);
        const auto code = fgmres_solve(singular_kernel, singular_diagonal, x.memptr(), rhs.memptr(), 1E-12);

        // oneMKL detects the zero norm of the generated vector and stops cleanly.
        REQUIRE(0 == code);
        REQUIRE(x.is_finite());
        REQUIRE_THAT(arma::norm(x, 2), Catch::Matchers::WithinAbs(0., 1E-8));
    }

    SECTION("solves repeatedly without leaking state between calls") {
        for(uword i = 0; i < 4; ++i) {
            vec x(n, arma::fill::zeros);
            MKL_INT iterations = -1;
            const vec scaled = static_cast<double>(i + 1) * expected;
            vec rhs = a * scaled;

            REQUIRE(0 == fgmres_solve(kernel, diagonal, x.memptr(), rhs.memptr(), 1E-12, &iterations));
            REQUIRE(iterations > 0);
            REQUIRE_THAT(max_difference(x, scaled), Catch::Matchers::WithinAbs(0., 1E-8));
        }
    }
}

//! @}
