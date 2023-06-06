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

#ifdef SUANPAN_MKL

#include <cstdio>
#include <mpi.h>
#include <memory>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int NUM_NODE = 4;

    MPI_Comm comm;
    MPI_Comm_spawn("solver.pardiso", MPI_ARGV_NULL, NUM_NODE, MPI_INFO_NULL, 0, MPI_COMM_SELF, &comm, MPI_ERRCODES_IGNORE);

    int iparm[64 + 7] = {0};
    int config[7] = {0};

    config[0] = 11;
    config[1] = 1;
    config[2] = 1;
    config[3] = 1;
    config[4] = 0;
    config[5] = 5;
    config[6] = 13;

    const auto n = config[5];
    const auto nnz = config[6];

    std::unique_ptr<int[]> ia(new int[n + 1]);
    std::unique_ptr<int[]> ja(new int[nnz]);
    std::unique_ptr<double[]> a(new double[nnz]);
    std::unique_ptr<double[]> b(new double[n]);

    ia[0] = 1;
    ia[1] = 4;
    ia[2] = 6;
    ia[3] = 9;
    ia[4] = 12;
    ia[5] = 14;
    ja[0] = 1;
    ja[1] = 2;
    ja[2] = 4;
    ja[3] = 1;
    ja[4] = 2;
    ja[5] = 3;
    ja[6] = 4;
    ja[7] = 5;
    ja[8] = 1;
    ja[9] = 3;
    ja[10] = 4;
    ja[11] = 2;
    ja[12] = 5;
    a[0] = 1.0;
    a[1] = -1.0;
    a[2] = -3.0;
    a[3] = -2.0;
    a[4] = 5.0;
    a[5] = 4.0;
    a[6] = 6.0;
    a[7] = 4.0;
    a[8] = -4.0;
    a[9] = 2.0;
    a[10] = 7.0;
    a[11] = 8.0;
    a[12] = -5.0;

    for(int i = 0; i < n; i++) b[i] = 1.0;

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

    std::unique_ptr<MPI_Request[]> requests(new MPI_Request[NUM_NODE + 5]);
    for(auto I = 0; I < NUM_NODE; ++I) MPI_Isend(&config, 7, MPI_INT, I, 0, comm, &requests[I]);
    MPI_Isend(&iparm, 64, MPI_INT, 0, 0, comm, &requests[NUM_NODE]);
    MPI_Isend(ia.get(), n + 1, MPI_INT, 0, 0, comm, &requests[NUM_NODE + 1]);
    MPI_Isend(ja.get(), nnz, MPI_INT, 0, 0, comm, &requests[NUM_NODE + 2]);
    MPI_Isend(a.get(), nnz, MPI_DOUBLE, 0, 0, comm, &requests[NUM_NODE + 3]);
    MPI_Isend(b.get(), n, MPI_DOUBLE, 0, 0, comm, &requests[NUM_NODE + 4]);

    MPI_Waitall(NUM_NODE + 5, requests.get(), MPI_STATUSES_IGNORE);

    MPI_Recv(b.get(), n, MPI_DOUBLE, 0, 0, comm, MPI_STATUS_IGNORE);

    MPI_Finalize();

    // for(int i = 0; i < n; i++) printf("x[%d] = %f\n", i, b[i]);

    return 0;
}

#else

int main(int argc, char* argv[]) { return 0; }

#endif
