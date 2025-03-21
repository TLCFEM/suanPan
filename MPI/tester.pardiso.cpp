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

#ifdef SUANPAN_MKL

#include <mpl/mpl.hpp>

void run() {
    constexpr int NUM_NODE = 2;

    const auto& comm_world{mpl::environment::comm_world()};
    const auto worker = comm_world.spawn(0, NUM_NODE, {"solver.pardiso"});
    const auto all = mpl::communicator(worker, mpl::communicator::order_low);

    int iparm[64]{};
    int config[7]{};

    iparm[0] = 1;   /* Solver default parameters overriden with provided by iparm */
    iparm[1] = 3;   /* Use METIS for fill-in reordering */
    iparm[5] = 0;   /* Write solution into x */
    iparm[7] = 2;   /* Max number of iterative refinement steps */
    iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
    iparm[12] = 1;  /* Switch on Maximum Weighted Matching algorithm (default for non-symmetric) */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[26] = 0;  /* Check input data for correctness */
    iparm[39] = 0;  /* Input: matrix/rhs/solution stored on master */

    config[0] = 11; // mtype
    config[1] = 1;  // nrhs
    config[2] = 1;  // maxfct
    config[3] = 1;  // mnum
    config[4] = 0;  // msglvl
    config[5] = 5;  // n
    config[6] = 13; // nnz

    const auto n = config[5];
    const auto nnz = config[6];

    std::vector<int> ia(n + 1), ja(nnz);
    std::vector<double> a(nnz), b(n, 1.);

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

    all.bcast(0, config);
    all.bcast(0, iparm);

    mpl::irequest_pool requests;

    requests.push(worker.isend(ia, 0, mpl::tag_t{0}));
    requests.push(worker.isend(ja, 0, mpl::tag_t{1}));
    requests.push(worker.isend(a, 0, mpl::tag_t{2}));
    requests.push(worker.isend(b, 0, mpl::tag_t{3}));

    requests.waitall();

    int error = -1;
    worker.recv(error, 0);
    if(0 == error) worker.recv(b, 0);

    for(int i = 0; i < n; i++) printf("x[%d] = %+.6f\n", i, b[i]);
}

int main(int argc, char* argv[]) {
    for(auto I = 0; I < 4; ++I) {
        printf("run %d\n", I);
        run();
    }

    return 0;
}

#else

int main(int, char**) { return 0; }

#endif
