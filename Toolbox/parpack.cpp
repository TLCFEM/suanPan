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
#include "parpack.h"

#include <Domain/MetaMat/operator_times.hpp>
#include <mpi.h>

using mat_ptr = std::shared_ptr<MetaMat<double>>;

#ifdef __cplusplus
extern "C" {
#endif
void pdsaupd_(const MPI_Fint*, blas_int*, char*, blas_int*, const char*, blas_int*, double*, double*, blas_int*, double*, blas_int*, blas_int*, blas_int*, double*, double*, blas_int*, blas_int*);
void pdseupd_(const MPI_Fint*, blas_int*, const char*, blas_int*, double*, double*, blas_int*, double*, const char*, blas_int*, const char*, blas_int*, double*, double*, blas_int*, double*, blas_int*, blas_int*, blas_int*, double*, double*, blas_int*, blas_int*);
void pdnaupd_(const MPI_Fint*, blas_int*, char*, blas_int*, const char*, blas_int*, double*, double*, blas_int*, double*, blas_int*, blas_int*, blas_int*, double*, double*, blas_int*, blas_int*);
void pdneupd_(const MPI_Fint*, blas_int*, const char*, blas_int*, double*, double*, double*, blas_int*, double*, double*, double*, const char*, blas_int*, const char*, blas_int*, double*, double*, blas_int*, double*, blas_int*, blas_int*, blas_int*, double*, double*, blas_int*, blas_int*);
#ifdef __cplusplus
}
#endif

namespace {
    // PARPACK stores a block of each vector on every rank, so its routines are given the local
    // dimension. Global quantities are recovered with the collectives below.
    class block_distribution {
        static constexpr int root{0};

        blas_int global_size{0};
        blas_int local_size{0};
        blas_int local_offset{0};
        std::vector<int> counts, displacements;

    public:
        explicit block_distribution(const blas_int N)
            : global_size(N) {
            int size, rank;
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            const auto base = N / size;
            const auto remainder = N % size;

            local_size = base + (rank < remainder ? 1 : 0);
            local_offset = rank * base + (rank < remainder ? rank : remainder);

            counts.resize(size);
            displacements.resize(size);

            int displacement = 0;
            for(auto I{0}; I < size; ++I) {
                counts[I] = static_cast<int>(base + (I < remainder ? 1 : 0));
                displacements[I] = displacement;
                displacement += counts[I];
            }
        }

        [[nodiscard]] blas_int global() const { return global_size; }
        [[nodiscard]] blas_int local() const { return local_size; }
        [[nodiscard]] blas_int offset() const { return local_offset; }

        [[nodiscard]] int minimum_local_size() const {
            int min_size = static_cast<int>(local_size);
            MPI_Allreduce(MPI_IN_PLACE, &min_size, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            return min_size;
        }

        void allgather(const double* local_ptr, double* global_ptr) const {
            MPI_Allgatherv(local_ptr, static_cast<int>(local_size), MPI_DOUBLE, global_ptr, counts.data(), displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        }

        void scatter(const double* global_ptr, double* local_ptr) const {
            MPI_Scatterv(global_ptr, counts.data(), displacements.data(), MPI_DOUBLE, local_ptr, static_cast<int>(local_size), MPI_DOUBLE, root, MPI_COMM_WORLD);
        }
    };

    vec collect_vector(const block_distribution& distribution, const double* local_ptr) {
        vec global(distribution.global());
        distribution.allgather(local_ptr, global.memptr());
        return global;
    }

    void distribute_vector(const block_distribution& distribution, const vec& global, double* local_ptr) { distribution.scatter(global.memptr(), local_ptr); }
} // namespace

int eig_solve(vec& eigval, mat& eigvec, const mat_ptr& K, const mat_ptr& M, const unsigned num, const char* WHICH) {
    static auto BMAT{'G'}; // generalized eigenvalue problem A*x=lambda*B*x

    const block_distribution distribution(K->n_cols);

    const auto COMM = MPI_Comm_c2f(MPI_COMM_WORLD);

    blas_int IDO{0}, INFO{0};
    const auto N = distribution.global();
    auto NLOCAL = distribution.local();
    auto NEV = std::min(static_cast<blas_int>(num), N - 1);
    auto TOL{0.};
    auto NCV = std::min(3 * NEV, static_cast<blas_int>(distribution.minimum_local_size())); // use a larger NCV to ensure convergence
    auto LWORKL = NCV * (NCV + 8);

    if(NEV >= NCV || NLOCAL <= 0) {
        suanpan_error("The problem is too small for the number of MPI processes.\n");
        return SUANPAN_FAIL;
    }

    blas_int IPARAM[11]{}, IPNTR[14]{};
    podarray<double> RESID(NLOCAL), V(uword(NLOCAL) * uword(NCV)), WORKD(3u * uword(NLOCAL)), WORKL(LWORKL);

    IPARAM[0] = 1;    // exact shift
    IPARAM[2] = 1000; // maximum iteration
    IPARAM[6] = 4;    // mode 4: K*x=lambda*KG*x

    // we choose shift to be -1 to ensure (M+K) is invertible
    // since we know the eigenvalues are all positive
    M += K;

    while(99 != IDO) {
        pdsaupd_(&COMM, &IDO, &BMAT, &NLOCAL, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &NLOCAL, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
        if(0 != INFO) break;
        if(-1 == IDO) {
            auto X = collect_vector(distribution, WORKD.memptr() + IPNTR[0] - 1);
            X = K * X;
            vec Y(N);
            INFO = M->solve(Y, X);
            if(0 != INFO) break;
            distribute_vector(distribution, Y, WORKD.memptr() + IPNTR[1] - 1);
        }
        else if(1 == IDO) {
            auto X = collect_vector(distribution, WORKD.memptr() + IPNTR[2] - 1);
            vec Y(N);
            INFO = M->solve(Y, X);
            if(0 != INFO) break;
            distribute_vector(distribution, Y, WORKD.memptr() + IPNTR[1] - 1);
        }
        else if(2 == IDO) {
            auto X = collect_vector(distribution, WORKD.memptr() + IPNTR[0] - 1);
            const auto Y = K * X;
            distribute_vector(distribution, Y, WORKD.memptr() + IPNTR[1] - 1);
        }
    }

    if(0 != INFO) {
        suanpan_error("Error code {} received.\n", INFO);
        return SUANPAN_FAIL;
    }

    suanpan_debug("Arnoldi iteration counter: {}.\n", IPARAM[2]);

    static blas_int RVEC{1};
    static auto HOWMNY{'A'};
    static auto SIGMA{-1.};

    podarray<blas_int> SELECT(NCV);
    podarray<double> Z(uword(NLOCAL) * uword(NEV));

    eigval.set_size(NEV);
    eigvec.set_size(static_cast<uword>(N), static_cast<uword>(NEV));

    pdseupd_(&COMM, &RVEC, &HOWMNY, SELECT.memptr(), eigval.memptr(), Z.memptr(), &NLOCAL, &SIGMA, &BMAT, &NLOCAL, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &NLOCAL, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

    if(0 == INFO)
        for(uword I{0}; I < uword(NEV); ++I) distribution.allgather(Z.memptr() + uword(NLOCAL) * I, eigvec.memptr() + uword(N) * I);

    return INFO;
}

int eig_solve(cx_vec& eigval, cx_mat& eigvec, const mat_ptr& K, const mat_ptr& M, const unsigned num, const char* WHICH) {
    static auto BMAT{'G'}; // generalized eigenvalue problem A*x=lambda*B*x

    const block_distribution distribution(K->n_rows);

    const auto COMM = MPI_Comm_c2f(MPI_COMM_WORLD);

    blas_int IDO{0}, INFO{0};
    const auto N = distribution.global();
    auto NLOCAL = distribution.local();
    auto NEV = std::min(static_cast<blas_int>(num), N - 2);
    auto TOL{0.};
    auto NCV = std::min(std::max(NEV + 3, 3 * NEV + 1), static_cast<blas_int>(distribution.minimum_local_size())); // use a larger NCV to ensure convergence
    auto LWORKL = 3 * NCV * (NCV + 2);

    if(NEV + 1 >= NCV || NLOCAL <= 0) {
        suanpan_error("The problem is too small for the number of MPI processes.\n");
        return SUANPAN_FAIL;
    }

    blas_int IPARAM[11]{}, IPNTR[14]{};
    podarray<double> RESID(NLOCAL), V(uword(NLOCAL) * uword(NCV)), WORKD(3u * uword(NLOCAL)), WORKL(LWORKL);

    IPARAM[0] = 1;    // exact shift
    IPARAM[2] = 1000; // maximum iteration
    IPARAM[6] = 3;    // mode 3: K*x=lambda*M*x

    while(99 != IDO) {
        pdnaupd_(&COMM, &IDO, &BMAT, &NLOCAL, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &NLOCAL, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
        if(0 != INFO) break;
        if(-1 == IDO) {
            auto X = collect_vector(distribution, WORKD.memptr() + IPNTR[0] - 1);
            X = M * X;
            vec Y(N);
            INFO = K->solve(Y, X);
            if(0 != INFO) break;
            distribute_vector(distribution, Y, WORKD.memptr() + IPNTR[1] - 1);
        }
        else if(1 == IDO) {
            auto X = collect_vector(distribution, WORKD.memptr() + IPNTR[2] - 1);
            vec Y(N);
            INFO = K->solve(Y, X);
            if(0 != INFO) break;
            distribute_vector(distribution, Y, WORKD.memptr() + IPNTR[1] - 1);
        }
        else if(2 == IDO) {
            auto X = collect_vector(distribution, WORKD.memptr() + IPNTR[0] - 1);
            const auto Y = M * X;
            distribute_vector(distribution, Y, WORKD.memptr() + IPNTR[1] - 1);
        }
    }

    if(0 != INFO) {
        suanpan_error("Error code {} received.\n", INFO);
        return SUANPAN_FAIL;
    }

    suanpan_debug("Arnoldi iteration counter: {}.\n", IPARAM[2]);

    static blas_int RVEC{1};
    static auto HOWMNY{'A'};
    static auto SIGMA{0.};

    podarray<blas_int> SELECT(NCV);
    podarray<double> DR(NEV + 1u), DI(NEV + 1u), Z(uword(NLOCAL) * uword(NEV + 1u)), WORKEV(3u * uword(NCV));

    pdneupd_(&COMM, &RVEC, &HOWMNY, SELECT.memptr(), DR.memptr(), DI.memptr(), Z.memptr(), &NLOCAL, &SIGMA, &SIGMA, WORKEV.memptr(), &BMAT, &NLOCAL, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &NLOCAL, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

    eigval.set_size(NEV);
    eigvec.set_size(static_cast<uword>(N), static_cast<uword>(NEV));

    podarray<double> GLOBAL_Z(uword(N) * uword(NEV + 1u));
    if(0 == INFO)
        for(uword I{0}; I < uword(NEV); ++I) distribution.allgather(Z.memptr() + uword(NLOCAL) * I, GLOBAL_Z.memptr() + uword(N) * I);

    // get eigenvalues
    for(uword I{0}; I < uword(NEV); ++I) eigval(I) = std::complex(DR(I), DI(I));

    // get eigenvectors
    for(uword I{0}; I < uword(NEV); ++I) {
        if(std::cmp_less(I + 1, NEV) && eigval[I] == std::conj(eigval[I + 1u])) {
            for(uword J{0}; J < uword(N); ++J) {
                eigvec.at(J, I) = std::complex(GLOBAL_Z[uword(N) * I + J], GLOBAL_Z[uword(N) * I + uword(N) + J]);
                eigvec.at(J, I + 1u) = std::complex(GLOBAL_Z[uword(N) * I + J], -GLOBAL_Z[uword(N) * I + uword(N) + J]);
            }
            ++I;
        }
        else if(std::cmp_equal(I + 1, NEV) && std::complex(eigval[I]).imag() != 0.)
            for(auto J = 0; J < N; ++J) eigvec.at(J, I) = std::complex(GLOBAL_Z[uword(N) * I + J], GLOBAL_Z[uword(N) * I + uword(N) + J]);
        else
            for(auto J = 0; J < N; ++J) eigvec.at(J, I) = std::complex(GLOBAL_Z[uword(N) * I + J], 0.);
    }

    return INFO;
}

int eig_solve(cx_vec& eigval, const std::shared_ptr<MetaMat<double>>& K, const unsigned num) {
    static auto BMAT{'I'};              // standard eigenvalue problem A*x=lambda*x
    static constexpr auto WHICH = "SR"; // smallest real part

    const block_distribution distribution(K->n_rows);

    const auto COMM = MPI_Comm_c2f(MPI_COMM_WORLD);

    blas_int IDO{0}, INFO{0};
    const auto N = distribution.global();
    auto NLOCAL = distribution.local();
    auto NEV = std::min(static_cast<blas_int>(num), N - 2);
    double TOL = std::numeric_limits<float>::epsilon();
    auto NCV = std::min(std::max(NEV + 3, 3 * NEV + 1), static_cast<blas_int>(distribution.minimum_local_size()));
    auto LWORKL = 3 * NCV * (NCV + 2);

    if(NEV + 1 >= NCV || NLOCAL <= 0) {
        suanpan_error("The problem is too small for the number of MPI processes.\n");
        return SUANPAN_FAIL;
    }

    blas_int IPARAM[11]{}, IPNTR[14]{};
    podarray<double> RESID(NLOCAL), V(uword(NLOCAL) * uword(NCV)), WORKD(3u * uword(NLOCAL)), WORKL(LWORKL);

    IPARAM[0] = 1;                           // exact shift
    IPARAM[2] = std::min(blas_int{1000}, N); // maximum iteration
    IPARAM[6] = 1;                           // mode 1: A*x=lambda*x

    while(99 != IDO) {
        pdnaupd_(&COMM, &IDO, &BMAT, &NLOCAL, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &NLOCAL, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
        if(0 != INFO) break;
        if(-1 == IDO || 1 == IDO) {
            const auto X = collect_vector(distribution, WORKD.memptr() + IPNTR[0] - 1);
            const auto Y = K * X;
            distribute_vector(distribution, Y, WORKD.memptr() + IPNTR[1] - 1);
        }
        else if(2 == IDO) {
            const auto X = collect_vector(distribution, WORKD.memptr() + IPNTR[0] - 1);
            distribute_vector(distribution, X, WORKD.memptr() + IPNTR[1] - 1);
        }
    }

    if(0 != INFO) {
        suanpan_error("Error code {} received.\n", INFO);
        return SUANPAN_FAIL;
    }

    suanpan_debug("Arnoldi iteration counter: {}.\n", IPARAM[2]);

    static blas_int RVEC{0};
    static auto HOWMNY{'A'};
    static auto SIGMA{0.};

    podarray<blas_int> SELECT(NCV);
    podarray<double> DR(NEV + 1u), DI(NEV + 1u), WORKEV(3u * uword(NCV));

    pdneupd_(&COMM, &RVEC, &HOWMNY, SELECT.memptr(), DR.memptr(), DI.memptr(), nullptr, &NLOCAL, &SIGMA, &SIGMA, WORKEV.memptr(), &BMAT, &NLOCAL, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &NLOCAL, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

    eigval.set_size(NEV);

    // get eigenvalues
    for(uword I{0}; I < uword(NEV); ++I) eigval(I) = std::complex(DR(I), DI(I));

    return INFO;
}
