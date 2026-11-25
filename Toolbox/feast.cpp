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

// ReSharper disable CppFunctionalStyleCast
// ReSharper disable IdentifierTypo
#include "feast.h"

#include <Domain/Factory.hpp>
#include <algorithm>
#include <array>

#ifdef __cplusplus
extern "C" {
#endif
void new_feastinit_(la_it*);
void new_dfeast_sygv_(const char*, la_it*, double*, la_it*, double*, la_it*, la_it*, double*, la_it*, double*, double*, la_it*, double*, double*, la_it*, double*, la_it*);
void new_dfeast_sbgv_(const char*, la_it*, la_it*, double*, la_it*, la_it*, double*, la_it*, la_it*, double*, la_it*, double*, double*, la_it*, double*, double*, la_it*, double*, la_it*);
void new_dfeast_scsrgv_(const char*, la_it*, double*, la_it*, la_it*, double*, la_it*, la_it*, la_it*, double*, la_it*, double*, double*, la_it*, double*, double*, la_it*, double*, la_it*);
#ifdef __cplusplus
}
#endif

namespace {
    constexpr char UPLO = 'F';

    std::array<la_it, 64> initialize_parameters() {
        std::array<la_it, 64> fpm{};

        new_feastinit_(fpm.data());

#ifdef SUANPAN_DEBUG
        fpm[0] = 1;
#endif

        fpm[41] = 0; // fpm(42): use full double precision

        return fpm;
    }

    csr_form<double, la_it> to_csr(const feast_mat& matrix) {
        if(matrix->is_empty()) return {};
        return csr_form<double, la_it>(to_triplet_form<double, la_it>(matrix), SparseBase::ONE);
    }

    int assign_result(vec& eigval, mat& eigvec, const vec& E, const mat& X, const la_it mode) {
        eigval = E.head(static_cast<uword>(mode));
        eigvec = X.head_cols(static_cast<uword>(mode));

        return SUANPAN_SUCCESS;
    }

    // FEAST stores a full band with the main diagonal on row `band` of an array with `2 * band + 1` rows.
    std::pair<std::vector<double>, std::vector<double>> to_full_band(const feast_mat& K, const feast_mat& M, const unsigned band) {
        const auto N = K->n_rows;
        const auto lda = 2 * band + 1;

        std::vector<double> A(lda * N, 0.), B(lda * N, 0.);

        for(uword J{0}; J < N; ++J)
            for(auto I = J > band ? J - band : 0; I <= std::min(N - 1, J + band); ++I) {
                const auto offset = band + I - J + J * lda;
                A[offset] = (*K)(I, J);
                B[offset] = (*M)(I, J);
            }

        return {std::move(A), std::move(B)};
    }

    // A lower band stores the main diagonal on row 0 of an array with `band + 1` rows.
    std::pair<std::vector<double>, std::vector<double>> to_lower_band(const feast_mat& K, const feast_mat& M, const unsigned band) {
        const auto N = K->n_rows;
        const auto lda = band + 1;

        std::vector<double> A(lda * N, 0.), B(lda * N, 0.);

        for(uword J{0}; J < N; ++J)
            for(auto I = J; I <= std::min(N - 1, J + band); ++I) {
                const auto offset = I - J + J * lda;
                A[offset] = (*K)(I, J);
                B[offset] = (*M)(I, J);
            }

        return {std::move(A), std::move(B)};
    }
} // namespace

int feast_solve(vec& eigval, mat& eigvec, const feast_mat& K, const feast_mat& M, const StorageScheme scheme, const unsigned num, const double centre, const double radius, const unsigned low_band, const unsigned up_band) {
    if(K->n_rows != K->n_cols || M->n_rows != M->n_cols || K->n_rows != M->n_rows) {
        suanpan_error("The system matrices must be square and of the same size.\n");
        return SUANPAN_FAIL;
    }

    if(M->is_empty()) {
        suanpan_error("The global mass matrix is empty, check if the element type supports mass and density settings.\n");
        return SUANPAN_FAIL;
    }

    auto fpm = initialize_parameters();

    auto N = static_cast<la_it>(K->n_rows);
    auto M0 = static_cast<la_it>(std::min<uword>(num, K->n_rows));

    if(0 == M0) {
        suanpan_error("The number of requested eigenvalues must be positive.\n");
        return SUANPAN_FAIL;
    }

    la_it mode{0}, loop{0}, info{0};
    double epsout{0.};

    double Emin = centre - radius;
    double Emax = centre + radius;

    vec E(M0, fill::zeros);
    mat X(N, M0, fill::zeros);
    vec R(M0, fill::zeros);

    switch(scheme) {
    case StorageScheme::FULL: {
        const auto size = static_cast<std::size_t>(N) * static_cast<std::size_t>(N);
        std::vector A(K->memptr(), K->memptr() + size);
        std::vector B(M->memptr(), M->memptr() + size);

        new_dfeast_sygv_(&UPLO, &N, A.data(), &N, B.data(), &N, fpm.data(), &epsout, &loop, &Emin, &Emax, &M0, E.memptr(), X.memptr(), &mode, R.memptr(), &info);
        break;
    }
    case StorageScheme::BAND: {
        const auto band = std::max(low_band, up_band);

        auto KL = static_cast<la_it>(band);
        auto LD = static_cast<la_it>(2 * band + 1);

        // SPIKE already stores a symmetric band matrix in the compact layout FEAST expects.
        auto* spike_k = dynamic_cast<BandMatSpike<double>*>(K.get());
        auto* spike_m = dynamic_cast<BandMatSpike<double>*>(M.get());

        if(nullptr != spike_k && nullptr != spike_m && low_band == up_band)
            new_dfeast_sbgv_(&UPLO, &N, &KL, spike_k->memptr(), &LD, &KL, spike_m->memptr(), &LD, fpm.data(), &epsout, &loop, &Emin, &Emax, &M0, E.memptr(), X.memptr(), &mode, R.memptr(), &info);
        else {
            suanpan_warning("The banded matrix is not stored in the SPIKE layout, converting it requires extra memory and impairs performance. Consider using the SPIKE solver directly.\n");

            auto [A, B] = to_full_band(K, M, band);

            new_dfeast_sbgv_(&UPLO, &N, &KL, A.data(), &LD, &KL, B.data(), &LD, fpm.data(), &epsout, &loop, &Emin, &Emax, &M0, E.memptr(), X.memptr(), &mode, R.memptr(), &info);
        }
        break;
    }
    case StorageScheme::BANDSYMM: {
        const auto band = std::max(low_band, up_band);
        auto [A, B] = to_lower_band(K, M, band);

        constexpr char LOWER = 'L';

        auto KL = static_cast<la_it>(band);
        auto LD = static_cast<la_it>(band + 1);

        new_dfeast_sbgv_(&LOWER, &N, &KL, A.data(), &LD, &KL, B.data(), &LD, fpm.data(), &epsout, &loop, &Emin, &Emax, &M0, E.memptr(), X.memptr(), &mode, R.memptr(), &info);
        break;
    }
    case StorageScheme::SPARSE:
    case StorageScheme::SPARSESYMM: {
        auto k = to_csr(K);
        auto m = to_csr(M);

        new_dfeast_scsrgv_(&UPLO, &N, k.val_mem(), k.row_mem(), k.col_mem(), m.val_mem(), m.row_mem(), m.col_mem(), fpm.data(), &epsout, &loop, &Emin, &Emax, &M0, E.memptr(), X.memptr(), &mode, R.memptr(), &info);
        break;
    }
    default: {
        suanpan_error("The current matrix storage scheme is not supported by the FEAST solver.\n");
        return SUANPAN_FAIL;
    }
    }

    if(0 != info) {
        suanpan_error("Error code {} received.\n", info);
        return SUANPAN_FAIL;
    }

    return assign_result(eigval, eigvec, E, X, mode);
}
