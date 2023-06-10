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

// ReSharper disable IdentifierTypo
#ifdef SUANPAN_MKL

#include <mpi.h>
#include <mkl_cluster_sparse_solver.h>
#include <memory>

int main(int argc, char** argv) {
    int error = 0, rank = -1;
    MPI_Comm parent, remote;
    int config[8];

    const auto mtype = &config[0];
    const auto nrhs = &config[1];
    const auto maxfct = &config[2];
    const auto mnum = &config[3];
    const auto msglvl = &config[4];
    const auto n = &config[5];
    const auto nnz = &config[6];
    const auto float_type = &config[7];

    std::unique_ptr<int[]> ia, ja;
    std::unique_ptr<double[]> a, b, x;

    auto finalise = [&] {
        if(0 == rank) MPI_Send(&error, 1, MPI_INT, 0, 0, parent);

        if(0 != error) MPI_Finalize();
        else {
            if(0 == rank) MPI_Send(x.get(), *n * *nrhs, *float_type > 0 ? MPI_DOUBLE : MPI_FLOAT, 0, 0, parent);
            error = MPI_Finalize();
        }

        return error;
    };

    error = MPI_Init(&argc, &argv);
    if(MPI_SUCCESS != error) return finalise();

    error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(MPI_SUCCESS != error) return finalise();

    error = MPI_Comm_get_parent(&parent);
    if(MPI_SUCCESS != error) return finalise();

    error = MPI_Intercomm_merge(parent, 0, &remote);
    if(MPI_SUCCESS != error) return finalise();

    error = MPI_Bcast(&config, 8, MPI_INT, 0, remote);
    if(MPI_SUCCESS != error) return finalise();

    const auto np = *n + 1;
    const auto nb = *n * *nrhs;

    int iparm[64];

    ia = std::make_unique<int[]>(np);
    ja = std::make_unique<int[]>(*nnz);
    a = std::make_unique<double[]>(*nnz);
    b = std::make_unique<double[]>(nb);
    x = std::make_unique<double[]>(nb);

    if(0 == rank) {
        MPI_Request requests[5];
        MPI_Irecv(&iparm, 64, MPI_INT, 0, 0, parent, &requests[0]);
        MPI_Irecv(ia.get(), np, MPI_INT, 0, 0, parent, &requests[1]);
        MPI_Irecv(ja.get(), *nnz, MPI_INT, 0, 0, parent, &requests[2]);
        MPI_Irecv(a.get(), *nnz, *float_type > 0 ? MPI_DOUBLE : MPI_FLOAT, 0, 0, parent, &requests[3]);
        MPI_Irecv(b.get(), nb, *float_type > 0 ? MPI_DOUBLE : MPI_FLOAT, 0, 0, parent, &requests[4]);
        MPI_Waitall(5, requests, MPI_STATUSES_IGNORE);

        iparm[0] = 1;                      /* Solver default parameters overriden with provided by iparm */
        iparm[1] = 2;                      /* Use METIS for fill-in reordering */
        iparm[5] = 0;                      /* Write solution into x */
        iparm[7] = 2;                      /* Max number of iterative refinement steps */
        iparm[9] = 13;                     /* Perturb the pivot elements with 1E-13 */
        iparm[10] = 1;                     /* Use nonsymmetric permutation and scaling MPS */
        iparm[12] = 1;                     /* Switch on Maximum Weighted Matching algorithm (default for non-symmetric) */
        iparm[17] = -1;                    /* Output: Number of nonzeros in the factor LU */
        iparm[18] = -1;                    /* Output: Mflops for LU factorization */
        iparm[26] = 0;                     /* Check input data for correctness */
        if(*float_type < 0) iparm[27] = 1; /* Single precision */
        iparm[39] = 0;                     /* Input: matrix/rhs/solution stored on master */
    }

    int64_t pt[64] = {0};

    // ReSharper disable once CppVariableCanBeMadeConstexpr
    const int comm = MPI_Comm_c2f(MPI_COMM_WORLD);

    int phase = 13;
    cluster_sparse_solver(&pt, maxfct, mnum, mtype, &phase, n, a.get(), ia.get(), ja.get(), nullptr, nrhs, iparm, msglvl, b.get(), x.get(), &comm, &error);
    if(0 != error) return finalise();

    phase = -1;
    cluster_sparse_solver(&pt, maxfct, mnum, mtype, &phase, n, nullptr, ia.get(), ja.get(), nullptr, nrhs, iparm, msglvl, nullptr, nullptr, &comm, &error);

    return finalise();
}

#else

int main(int argc, char** argv) { return 0; }

#endif
