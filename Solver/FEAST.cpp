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

#include "FEAST.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <feast/feast.h>

char FEAST::UPLO = 'F';

int FEAST::linear_solve(const shared_ptr<LongFactory>& W) const {
    auto& mass = W->get_mass();
    auto& stiffness = W->get_stiffness();

    std::vector fpm(64, 0);

    new_feastinit_(fpm.data());

#ifdef SUANPAN_DEBUG
    fpm[0] = 1;
#endif

    std::vector output(4, 0);
    std::vector input(4, 0.);
    input[1] = centre - radius; // centre
    input[2] = centre + radius; // radius

    output[1] = static_cast<int>(eigen_num);

    int N = static_cast<int>(W->get_size());

    auto M = static_cast<int>(eigen_num);
    std::vector R(M, 0.);
    std::vector E(M, 0.);
    M *= N;
    std::vector X(M, 0.);

    if(const auto scheme = W->get_storage_scheme(); StorageScheme::FULL == scheme) new_dfeast_sygv_(&UPLO, &N, stiffness->memptr(), &N, mass->memptr(), &N, fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
    else if(StorageScheme::SPARSE == scheme || StorageScheme::SPARSESYMM == scheme) {
        auto fs = std::async([&] {
            auto triplet_k = to_triplet_form<double, int>(stiffness);
            return csr_form<double, int>(triplet_k, SparseBase::ONE);
        });
        auto fm = std::async([&] {
            auto triplet_m = to_triplet_form<double, int>(mass);
            return csr_form<double, int>(triplet_m, SparseBase::ONE);
        });

        auto t_stiff = fs.get();
        auto t_mass = fm.get();

        new_dfeast_scsrgv_(&UPLO, &N, t_stiff.val_mem(), t_stiff.row_mem(), t_stiff.col_mem(), t_mass.val_mem(), t_mass.row_mem(), t_mass.col_mem(), fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
    }
    else if(StorageScheme::BAND == scheme || StorageScheme::BANDSYMM == scheme) {
        fpm[41] = 0;

        unsigned l, u;
        W->get_bandwidth(l, u);
        auto KL = static_cast<int>(l);
        const auto KU = static_cast<int>(u);
        auto LD = KL + KU + 1;

        new_dfeast_sbgv_(&UPLO, &N, &KL, stiffness->memptr(), &LD, &KL, mass->memptr(), &LD, fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
    }
    else {
        suanpan_error("The current matrix storage scheme is not supported by the FEAST solver.\n");
        return SUANPAN_FAIL;
    }

    if(0 != output[3]) {
        suanpan_error("Error code {} received.\n", output[3]);
        return SUANPAN_FAIL;
    }

    W->modify_eigenvalue() = vec(E.data(), output[2]);
    W->modify_eigenvector() = mat(X.data(), N, output[2]);

    return SUANPAN_SUCCESS;
}

int FEAST::quadratic_solve(const shared_ptr<LongFactory>& W) const {
    std::vector fpm(64, 0);

    new_feastinit_(fpm.data());

#ifdef SUANPAN_DEBUG
    fpm[0] = 1;
#endif

    int N = static_cast<int>(W->get_size());

    std::vector output(4, 0);
    std::vector input(4, 0.);
    input[0] = centre; // centre
    input[1] = 0.;     // centre
    input[2] = radius; // radius

    output[1] = static_cast<int>(eigen_num);

    auto M = 2 * static_cast<int>(eigen_num);
    std::vector R(M, 0.);
    std::vector E(M, 0.);
    M *= 2 * N;
    std::vector X(M, 0.);

    int P = 2;

    auto fk = std::async([&] {
        auto triplet_k = to_triplet_form<double, int>(W->get_stiffness());
        return csr_form<double, int>(triplet_k, SparseBase::ONE);
    });
    auto fd = std::async([&] {
        auto triplet_d = to_triplet_form<double, int>(W->get_damping());
        return csr_form<double, int>(triplet_d, SparseBase::ONE);
    });
    auto fm = std::async([&] {
        auto triplet_m = to_triplet_form<double, int>(W->get_mass());
        return csr_form<double, int>(triplet_m, SparseBase::ONE);
    });

    const auto t_stiff = fk.get();
    const auto t_damping = fd.get();
    const auto t_mass = fm.get();

    size_t n_elem1 = std::max(std::max(t_stiff.n_elem, t_damping.n_elem), t_mass.n_elem);
    size_t n_elem2 = n_elem1 + n_elem1;
    size_t n_elem3 = n_elem1 + n_elem2;
    size_t n_size1 = N;
    size_t n_size2 = n_size1 + n_size1;
    size_t n_size3 = n_size2 + n_size1;

    std::vector A(n_elem3, 0.);
    std::vector JA(n_elem3, 0);
    std::vector IA(n_size3, 0);

    suanpan::for_each(0, t_stiff.n_elem, [&](const int I) {
        A[I] = t_stiff.val_mem()[I];
        JA[I] = t_stiff.col_mem()[I];
    });
    suanpan::for_each(0, t_damping.n_elem, [&](const int I) {
        A[I + n_elem1] = t_damping.val_mem()[I];
        JA[I + n_elem1] = t_damping.col_mem()[I];
    });
    suanpan::for_each(0, t_mass.n_elem, [&](const int I) {
        A[I + n_elem2] = t_mass.val_mem()[I];
        JA[I + n_elem2] = t_mass.col_mem()[I];
    });

    suanpan::for_each(0, N, [&](const int I) {
        JA[I] = t_stiff.row_mem()[I];
        JA[I + n_size1] = t_damping.row_mem()[I];
        JA[I + n_size2] = t_mass.row_mem()[I];
    });

    new_dfeast_gcsrpev_(&P, &N, A.data(), IA.data(), JA.data(), fpm.data(), &input[3], &output[0], &input[0], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);

    if(0 != output[3]) {
        suanpan_error("Error code {} received.\n", output[3]);
        return SUANPAN_FAIL;
    }

    auto& eigval = W->modify_eigenvalue();
    eigval.set_size(output[2]);

    for(uword I = 0; I < eigval.n_elem; ++I) eigval(I) = E[2 * I];

    auto& eigvec = W->modify_eigenvector();
    eigvec.resize(N, output[2]);

    for(uword I = 0; I < eigvec.n_elem; ++I) eigvec(I) = X[2 * I];

    return SUANPAN_SUCCESS;
}

FEAST::FEAST(const unsigned T, const unsigned N, const double C, const double R, const bool Q)
    : Solver(T)
    , quadratic(Q)
    , eigen_num(N)
    , centre(C)
    , radius(R) {}

int FEAST::initialize() {
    auto& G = get_integrator();

    if(nullptr == G) {
        suanpan_error("A valid integrator is required.\n");
        return SUANPAN_FAIL;
    }

    auto& W = G->get_domain()->get_factory();

    if(const auto scheme = W->get_storage_scheme(); StorageScheme::SYMMPACK == scheme) {
        suanpan_error("The symmetric pack storage is not supported.\n");

        return SUANPAN_FAIL;
    }
    else if((StorageScheme::BAND == scheme || StorageScheme::BANDSYMM == scheme) && SolverType::SPIKE != W->get_solver_type()) {
        suanpan_error("The SPIKE system solver is required for banded storage.\n");

        return SUANPAN_FAIL;
    }

    return SUANPAN_SUCCESS;
}

int FEAST::analyze() {
    auto& G = get_integrator();
    const auto D = G->get_domain();
    auto& W = D->get_factory();

    if(SUANPAN_SUCCESS != G->process_modifier()) return SUANPAN_FAIL;

    D->assemble_trial_mass();
    D->assemble_trial_stiffness();
    if(quadratic) D->assemble_trial_damping();

    // if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;

    return quadratic ? quadratic_solve(W) : linear_solve(W);
}

void FEAST::print() {
    suanpan_info("An eigen solver using FEAST solver.\n");
}
