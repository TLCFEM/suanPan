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
#include "pfeast.h"

#include <Domain/Factory.hpp>
#include <algorithm>
#include <array>

#ifdef __cplusplus
extern "C" {
#endif
void new_pfeastinit_(la_it*, la_it*, la_it*);
void new_pdfeast_scsrgv_(const char*, la_it*, double*, la_it*, la_it*, double*, la_it*, la_it*, la_it*, double*, la_it*, double*, double*, la_it*, double*, double*, la_it*, double*, la_it*);
#ifdef __cplusplus
}
#endif

namespace {
    constexpr char UPLO = 'F';

    std::array<la_it, 64> initialize_parameters(const la_it size) {
        std::array<la_it, 64> fpm{};

        // Keep the L3 group no larger than the problem so that every rank owns at least one row.
        auto n_l3 = std::max<la_it>(1, std::min(static_cast<la_it>(comm_size), size));
        while(n_l3 > 1 && 0 != comm_size % n_l3) --n_l3;

        // Take the processor grid from the existing communicator.
        auto comm = static_cast<la_it>(MPI_Comm_c2f(comm_world.native_handle()));
        new_pfeastinit_(fpm.data(), &comm, &n_l3);

#ifdef SUANPAN_DEBUG
        fpm[0] = 1;
#endif

        fpm[41] = 0; // fpm(42): use full double precision
        fpm[42] = 1; // fpm(43): use the PIFEAST MPI drivers

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
} // namespace

int feast_solve(vec& eigval, mat& eigvec, const feast_mat& K, const feast_mat& M, const StorageScheme scheme, const unsigned num, const double centre, const double radius, const unsigned, const unsigned) {
    if(StorageScheme::SPARSE != scheme && StorageScheme::SPARSESYMM != scheme) {
        suanpan_error("Only sparse storage is supported by the parallel FEAST solver.\n");
        return SUANPAN_FAIL;
    }

    if(M->is_empty()) {
        suanpan_error("The global mass matrix is empty, check if the element type supports mass and density settings.\n");
        return SUANPAN_FAIL;
    }

    if(K->n_rows != K->n_cols || M->n_rows != M->n_cols || K->n_rows != M->n_rows) {
        suanpan_error("The system matrices must be square and of the same size.\n");
        return SUANPAN_FAIL;
    }

    auto N = static_cast<la_it>(K->n_rows);
    auto M0 = static_cast<la_it>(std::min<uword>(num, K->n_rows));

    if(0 == M0) {
        suanpan_error("The number of requested eigenvalues must be positive.\n");
        return SUANPAN_FAIL;
    }

    auto fpm = initialize_parameters(N);

    la_it mode{0}, loop{0}, info{0};
    double epsout{0.};

    double Emin = centre - radius;
    double Emax = centre + radius;

    vec E(M0, fill::zeros);
    mat X(N, M0, fill::zeros);
    vec R(M0, fill::zeros);

    // The assembly already provides the global matrices on every rank, PFEAST distributes the rows.
    auto k = to_csr(K);
    auto m = to_csr(M);

    new_pdfeast_scsrgv_(&UPLO, &N, k.val_mem(), k.row_mem(), k.col_mem(), m.val_mem(), m.row_mem(), m.col_mem(), fpm.data(), &epsout, &loop, &Emin, &Emax, &M0, E.memptr(), X.memptr(), &mode, R.memptr(), &info);

    if(0 != info) {
        suanpan_error("Error code {} received.\n", info);
        return SUANPAN_FAIL;
    }

    return assign_result(eigval, eigvec, E, X, mode);
}
