/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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
#include "arpack.h"

#include <Domain/MetaMat/operator_times.hpp>

using mat_ptr = std::shared_ptr<MetaMat<double>>;

int eig_solve(vec& eigval, mat& eigvec, const mat_ptr& K, const mat_ptr& M, const unsigned num, const char* WHICH) {
    static auto BMAT = 'G'; // generalized eigenvalue problem A*x=lambda*B*x

    blas_int IDO = 0, INFO = 0;
    auto N = static_cast<blas_int>(K->n_cols);
    auto NEV = std::min(static_cast<blas_int>(num), N - 1);
    auto TOL = 0.;
    auto NCV = std::min(3 * NEV, N); // use a larger NCV to ensure convergence
    auto LWORKL = 2 * NCV * (NCV + 8);

    blas_int IPARAM[11]{}, IPNTR[14]{};
    podarray<double> RESID(N), V(uword(N) * uword(NCV)), WORKD(5 * uword(N)), WORKL(LWORKL);

    IPARAM[0] = 1;    // exact shift
    IPARAM[2] = 1000; // maximum iteration
    IPARAM[6] = 4;    // mode 4: K*x=lambda*KG*x

    // we choose shift to be -1 to ensure (M+K) is invertible
    // since we know the eigenvalues are all positive
    M += K;

    while(99 != IDO) {
        arma_fortran(arma_dsaupd)(&IDO, &BMAT, &N, (char*)WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &N, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
        if(0 != INFO) break;
        // ReSharper disable once CppEntityAssignedButNoRead
        if(vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false, true); -1 == IDO) {
            vec X(WORKD.memptr() + IPNTR[0] - 1, N, false, true);
            X = K * X;
            INFO = M->solve(Y, X);
            if(0 != INFO) break;
        }
        else if(1 == IDO) {
            const vec X(WORKD.memptr() + IPNTR[2] - 1, N, false, true);
            INFO = M->solve(Y, X);
            if(0 != INFO) break;
        }
        else if(2 == IDO) {
            const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false, true);
            // ReSharper disable once CppDFAUnusedValue
            Y = K * X;
        }
    }

    if(0 != INFO) {
        suanpan_error("Error code {} received.\n", INFO);
        return SUANPAN_FAIL;
    }

    suanpan_debug("Arnoldi iteration counter: {}.\n", IPARAM[2]);

    static blas_int RVEC = 1;
    static auto HOWMNY = 'A';
    static auto SIGMA = -1.;

    podarray<blas_int> SELECT(NCV);

    eigval.set_size(NEV);
    eigvec.set_size(N, NEV);

    arma_fortran(arma_dseupd)(&RVEC, &HOWMNY, SELECT.memptr(), eigval.memptr(), eigvec.memptr(), &N, &SIGMA, &BMAT, &N, (char*)WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &N, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

    return INFO;
}

int eig_solve(cx_vec& eigval, cx_mat& eigvec, const mat_ptr& K, const mat_ptr& M, const unsigned num, const char* WHICH) {
    static auto BMAT = 'G'; // generalized eigenvalue problem A*x=lambda*B*x

    blas_int IDO = 0, INFO = 0;
    auto N = static_cast<blas_int>(K->n_rows);
    auto NEV = std::min(static_cast<blas_int>(num), N - 2);
    auto TOL = 0.;
    auto NCV = std::min(std::max(NEV + 3, 3 * NEV + 1), N); // use a larger NCV to ensure convergence
    auto LWORKL = 3 * NCV * (NCV + 2);

    blas_int IPARAM[11]{}, IPNTR[14]{};
    podarray<double> RESID(N), V(N * uword(NCV)), WORKD(3llu * N), WORKL(LWORKL);

    IPARAM[0] = 1;    // exact shift
    IPARAM[2] = 1000; // maximum iteration
    IPARAM[6] = 3;    // mode 3: K*x=lambda*M*x

    while(99 != IDO) {
        arma_fortran(arma_dnaupd)(&IDO, &BMAT, &N, (char*)WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &N, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
        if(0 != INFO) break;
        // ReSharper disable once CppEntityAssignedButNoRead
        if(vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false, true); -1 == IDO) {
            vec X(WORKD.memptr() + IPNTR[0] - 1, N, false, true);
            X = M * X;
            INFO = K->solve(Y, X);
            if(0 != INFO) break;
        }
        else if(1 == IDO) {
            const vec X(WORKD.memptr() + IPNTR[2] - 1, N, false, true);
            INFO = K->solve(Y, X);
            if(0 != INFO) break;
        }
        else if(2 == IDO) {
            const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false, true);
            // ReSharper disable once CppDFAUnusedValue
            Y = M * X;
        }
    }

    if(0 != INFO) {
        suanpan_error("Error code {} received.\n", INFO);
        return SUANPAN_FAIL;
    }

    suanpan_debug("Arnoldi iteration counter: {}.\n", IPARAM[2]);

    static blas_int RVEC = 1;
    static auto HOWMNY = 'A';
    static auto SIGMAR = 0., SIGMAI = 0.;

    podarray<blas_int> SELECT(NCV);
    podarray<double> DR(NEV + 1llu), DI(NEV + 1llu), Z(N * (NEV + 1llu)), WORKEV(3llu * NCV);

    arma_fortran(arma_dneupd)(&RVEC, &HOWMNY, SELECT.memptr(), DR.memptr(), DI.memptr(), Z.memptr(), &N, &SIGMAR, &SIGMAI, WORKEV.memptr(), &BMAT, &N, (char*)WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &N, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

    eigval.set_size(NEV);
    eigvec.set_size(N, NEV);

    // get eigenvalues
    for(uword I = 0; I < uword(NEV); ++I) eigval(I) = std::complex(DR(I), DI(I));

    // get eigenvectors
    for(uword I = 0; I < uword(NEV); ++I) {
        if(I < NEV - 1llu && eigval[I] == std::conj(eigval[I + 1llu])) {
            for(uword J = 0; J < uword(N); ++J) {
                eigvec.at(J, I) = std::complex(Z[N * I + J], Z[N * I + N + J]);
                eigvec.at(J, I + 1llu) = std::complex(Z[N * I + J], -Z[N * I + N + J]);
            }
            ++I;
        }
        else if(I == NEV - 1llu && std::complex(eigval[I]).imag() != 0.)
            for(auto J = 0; J < N; ++J) eigvec.at(J, I) = std::complex(Z[N * I + J], Z[N * I + N + J]);
        else
            for(auto J = 0; J < N; ++J) eigvec.at(J, I) = std::complex(Z[N * I + J], 0.);
    }

    return INFO;
}

int eig_solve(cx_vec& eigval, const std::shared_ptr<MetaMat<double>>& K, const unsigned num) {
    static auto BMAT = 'I';             // standard eigenvalue problem A*x=lambda*x
    static constexpr auto WHICH = "SR"; // smallest real part

    blas_int IDO = 0, INFO = 0;
    auto N = static_cast<blas_int>(K->n_rows);
    auto NEV = std::min(static_cast<blas_int>(num), N - 2);
    double TOL = std::numeric_limits<float>::epsilon();
    auto NCV = std::min(std::max(NEV + 3, 3 * NEV + 1), N);
    auto LWORKL = 3 * NCV * (NCV + 2);

    blas_int IPARAM[11]{}, IPNTR[14]{};
    podarray<double> RESID(N), V(N * uword(NCV)), WORKD(3llu * N), WORKL(LWORKL);

    IPARAM[0] = 1;    // exact shift
    IPARAM[2] = 1000; // maximum iteration
    IPARAM[6] = 1;    // mode 1: A*x=lambda*x

    while(99 != IDO) {
        arma_fortran(arma_dnaupd)(&IDO, &BMAT, &N, (char*)WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &N, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
        if(0 != INFO) break;
        const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false, true);
        if(vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false, true); -1 == IDO || 1 == IDO) Y = K * X;
        else if(2 == IDO) Y = X;
    }

    if(0 != INFO) {
        suanpan_error("Error code {} received.\n", INFO);
        return SUANPAN_FAIL;
    }

    suanpan_debug("Arnoldi iteration counter: {}.\n", IPARAM[2]);

    static blas_int RVEC = 0;
    static auto HOWMNY = 'A';
    static auto SIGMAR = 0., SIGMAI = 0.;

    podarray<blas_int> SELECT(NCV);
    podarray<double> DR(NEV + 1llu), DI(NEV + 1llu), WORKEV(3llu * NCV);

    arma_fortran(arma_dneupd)(&RVEC, &HOWMNY, SELECT.memptr(), DR.memptr(), DI.memptr(), nullptr, &N, &SIGMAR, &SIGMAI, WORKEV.memptr(), &BMAT, &N, (char*)WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &N, IPARAM, IPNTR, WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

    eigval.set_size(NEV);

    // get eigenvalues
    for(uword I = 0; I < uword(NEV); ++I) eigval(I) = std::complex(DR(I), DI(I));

    return INFO;
}
