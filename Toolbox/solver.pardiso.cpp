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

int finalise(int64_t error) {
    if(0 == error) error = MPI_Finalize();
    else MPI_Finalize();
    return static_cast<int>(error);
}

int main(int argc, char** argv) {
    int error = MPI_Init(&argc, &argv);
    if(MPI_SUCCESS != error) return finalise(error);

    int rank;
    error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(MPI_SUCCESS != error) return finalise(error);

    MPI_Comm parent, remote;
    error = MPI_Comm_get_parent(&parent);
    if(MPI_SUCCESS != error) return finalise(error);
    error = MPI_Intercomm_merge(parent, 0, &remote);
    if(MPI_SUCCESS != error) return finalise(error);

    int config[7];
    error = MPI_Bcast(&config, 7, MPI_INT, 0, remote);
    if(MPI_SUCCESS != error) return finalise(error);

    const auto mtype = &config[0];
    const auto nrhs = &config[1];
    const auto maxfct = &config[2];
    const auto mnum = &config[3];
    const auto msglvl = &config[4];
    const auto n = &config[5];
    const auto nnz = config[6];
    const auto np = *n + 1;

    int iparm[64];

    const std::unique_ptr<int[]> ia(new int[np]), ja(new int[nnz]);
    const std::unique_ptr<double[]> a(new double[nnz]), b(new double[*n]), x(new double[*n]);

    if(0 == rank) {
        std::unique_ptr<MPI_Request[]> requests(new MPI_Request[5]);
        MPI_Irecv(&iparm, 64, MPI_INT, 0, 0, parent, &requests[0]);
        MPI_Irecv(ia.get(), np, MPI_INT, 0, 0, parent, &requests[1]);
        MPI_Irecv(ja.get(), nnz, MPI_INT, 0, 0, parent, &requests[2]);
        MPI_Irecv(a.get(), nnz, MPI_DOUBLE, 0, 0, parent, &requests[3]);
        MPI_Irecv(b.get(), *n, MPI_DOUBLE, 0, 0, parent, &requests[4]);
        MPI_Waitall(5, requests.get(), MPI_STATUSES_IGNORE);

        iparm[0] = 1;   /* Solver default parameters overriden with provided by iparm */
        iparm[1] = 2;   /* Use METIS for fill-in reordering */
        iparm[5] = 0;   /* Write solution into x */
        iparm[7] = 2;   /* Max number of iterative refinement steps */
        iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
        iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
        iparm[12] = 1;  /* Switch on Maximum Weighted Matching algorithm (default for non-symmetric) */
        iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
        iparm[18] = -1; /* Output: Mflops for LU factorization */
        iparm[26] = 0;  /* Check input data for correctness */
        iparm[39] = 0;  /* Input: matrix/rhs/solution stored on master */
    }

    int64_t pt[64] = {0};

    // ReSharper disable once CppVariableCanBeMadeConstexpr
    const int comm = MPI_Comm_c2f(MPI_COMM_WORLD);

    int phase = 13;
    cluster_sparse_solver(&pt, maxfct, mnum, mtype, &phase, n, a.get(), ia.get(), ja.get(), nullptr, nrhs, iparm, msglvl, b.get(), x.get(), &comm, &error);
    if(0 != error) return finalise(error);

    if(0 == rank) MPI_Send(x.get(), *n, MPI_DOUBLE, 0, 0, parent);

    phase = -1;
    cluster_sparse_solver(&pt, maxfct, mnum, mtype, &phase, n, nullptr, ia.get(), ja.get(), nullptr, nrhs, iparm, msglvl, nullptr, nullptr, &comm, &error);

    return finalise(error);
}

#else

int main(int argc, char** argv) { return 0; }

#endif
