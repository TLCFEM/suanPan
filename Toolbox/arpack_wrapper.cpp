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

// ReSharper disable CppFunctionalStyleCast
// ReSharper disable IdentifierTypo
#include "arpack_wrapper.h"
#include <Domain/MetaMat/operator_times.hpp>

int eig_solve(vec& eigval, mat& eigvec, const std::shared_ptr<MetaMat<double>>& K, const std::shared_ptr<MetaMat<double>>& M, const unsigned num, const char* form) {
    auto IDO = 0;
    auto BMAT = 'G'; // generalized eigenvalue problem A*x=lambda*M*x
    auto N = static_cast<int>(K->n_cols);
    char WHICH[2];
    for(auto I = 0; I < 2; ++I) WHICH[I] = form[I];
    auto NEV = std::min(static_cast<int>(num), N - 1);
    auto TOL = 0.;
    auto NCV = std::min(2 * NEV, N);
    auto LDV = N;
    auto LWORKL = 2 * NCV * (NCV + 8);
    auto INFO = 0;

    podarray<int> IPARAM(11), IPNTR(14);
    podarray<double> RESID(N), V(uword(N) * uword(NCV)), WORKD(5 * uword(N)), WORKL(LWORKL);

    IPARAM(0) = 1;    // exact shift
    IPARAM(2) = 1000; // maximum iteration
    IPARAM(6) = 4;    // mode 4: K*x=lambda*KG*x

    auto SIGMA = -1.;

    M += K;

    while(99 != IDO) {
        arma_fortran(arma_dsaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
        // ReSharper disable once CppEntityAssignedButNoRead
        if(vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false); -1 == IDO) {
            vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
            X = K * X;
            Y = M->solve(X);
        }
        else if(1 == IDO) {
            const vec X(WORKD.memptr() + IPNTR[2] - 1, N, false);
            Y = M->solve(X);
        }
        else if(2 == IDO) {
            const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
            Y = K * X;
        }
        else if(0 != INFO) break;
    }

    if(0 != INFO) {
        SP_E("Error code {} received.\n", INFO);
        return SUANPAN_FAIL;
    }

    SP_D("Arnoldi iteration counter: {}.\n", IPARAM(2));

    auto RVEC = 1;
    auto HOWMNY = 'A';
    auto LDZ = N;

    podarray<int> SELECT(NCV);

    eigval.set_size(NEV);
    eigvec.set_size(N, NEV);

    arma_fortran(arma_dseupd)(&RVEC, &HOWMNY, SELECT.memptr(), eigval.memptr(), eigvec.memptr(), &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

    return INFO;
}

int eig_solve(cx_vec& eigval, cx_mat& eigvec, const std::shared_ptr<MetaMat<double>>& K, const std::shared_ptr<MetaMat<double>>& M, const unsigned num, const char* form) {
    auto IDO = 0;
    auto BMAT = 'G'; // standard eigenvalue problem A*x=lambda*x
    auto N = static_cast<int>(K->n_rows);
    char WHICH[2];
    for(auto I = 0; I < 2; ++I) WHICH[I] = form[I];
    auto NEV = std::min(static_cast<int>(num), N - 2);
    auto TOL = 0.;
    auto NCV = std::min(std::max(NEV + 3, 2 * NEV + 1), N);
    auto LDV = N;
    auto LWORKL = 3 * NCV * (NCV + 2);
    auto INFO = 0;

    podarray<int> IPARAM(11), IPNTR(14);
    podarray<double> RESID(N), V(N * uword(NCV)), WORKD(3llu * N), WORKL(LWORKL);

    IPARAM(0) = 1;    // exact shift
    IPARAM(2) = 1000; // maximum iteration
    IPARAM(6) = 3;    // mode 3: K*x=lambda*KG*x

    auto SIGMAR = 0.;
    auto SIGMAI = 0.;

    // K += M;

    while(99 != IDO) {
        arma_fortran(arma_dnaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
        // ReSharper disable once CppEntityAssignedButNoRead
        if(vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false); -1 == IDO) {
            vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
            X = M * X;
            Y = K->solve(X);
        }
        else if(1 == IDO) {
            const vec X(WORKD.memptr() + IPNTR[2] - 1, N, false);
            Y = K->solve(X);
        }
        else if(2 == IDO) {
            const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
            Y = M * X;
        }
        else if(0 != INFO) break;
    }

    if(0 != INFO) {
        SP_E("Error code {} received.\n", INFO);
        return SUANPAN_FAIL;
    }

    auto RVEC = 1;
    auto HOWMNY = 'A';
    auto LDZ = N;

    podarray<int> SELECT(NCV);
    podarray<double> DR(NEV + 1llu), DI(NEV + 1llu), Z(N * (NEV + 1llu)), WORKEV(3llu * NCV);

    arma_fortran(arma_dneupd)(&RVEC, &HOWMNY, SELECT.memptr(), DR.memptr(), DI.memptr(), Z.memptr(), &LDZ, &SIGMAR, &SIGMAI, WORKEV.memptr(), &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

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
        else if(I == NEV - 1llu && std::complex(eigval[I]).imag() != 0.) for(auto J = 0; J < N; ++J) eigvec.at(J, I) = std::complex(Z[N * I + J], Z[N * I + N + J]);
        else for(auto J = 0; J < N; ++J) eigvec.at(J, I) = std::complex(Z[N * I + J], 0.);
    }

    return INFO;
}
