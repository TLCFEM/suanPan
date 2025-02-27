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
/**
 * @class SparseMatMPIPARDISO
 * @brief A SparseMatMPI class that holds matrices.
 *
 * * MUMPS uses int.
 *
 * @author tlc
 * @date 27/02/2025
 * @version 0.1.0
 * @file SparseMatMPI.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppClangTidyClangDiagnosticMissingFieldInitializers
#ifndef SPARSEMATMPI_HPP
#define SPARSEMATMPI_HPP

#if defined(SUANPAN_MPI) && defined(SUANPAN_MKL)

#include <mpl/mpl.hpp>

extern int SUANPAN_NUM_NODES;

template<sp_d T> class SparseMatMPIPARDISO final : public SparseMat<T> {
    int iparm[64]{};

protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatMPIPARDISO(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) {
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
        iparm[34] = 1;  /* zero-based indexing */
        iparm[39] = 0;  /* Input: matrix/rhs/solution stored on master */
    }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatMPIPARDISO>(*this); }
};

template<sp_d T> int SparseMatMPIPARDISO<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    X.set_size(B.n_rows, B.n_cols);

    csr_form<T, int> csr_mat(this->triplet_mat, SparseBase::ZERO, true);

    const auto n = static_cast<int>(B.n_rows);
    const auto nrhs = static_cast<int>(B.n_cols);
    const auto nnz = static_cast<int>(csr_mat.n_elem);

    const auto& comm_world{mpl::environment::comm_world()};
    const auto worker = comm_world.spawn(0, SUANPAN_NUM_NODES, {"solver.pardiso"});
    const auto all = mpl::communicator(worker, mpl::communicator::order_low);

    int config[7]{};

    config[0] = 11;   // mtype
    config[1] = nrhs; // nrhs
    config[2] = 1;    // maxfct
    config[3] = 1;    // mnum
    config[4] = 0;    // msglvl
    config[5] = n;    // n
    config[6] = nnz;  // nnz

    all.bcast(0, config);

    if(std::is_same_v<T, double>) iparm[27] = 0;
    else iparm[27] = 1;

    all.bcast(0, iparm);

    mpl::irequest_pool requests;

    requests.push(worker.isend(csr_mat.row_mem(), csr_mat.row_mem() + n + 1, 0, mpl::tag_t{0}));
    requests.push(worker.isend(csr_mat.col_mem(), csr_mat.col_mem() + nnz, 0, mpl::tag_t{1}));
    requests.push(worker.isend(csr_mat.val_mem(), csr_mat.val_mem() + nnz, 0, mpl::tag_t{2}));
    requests.push(worker.isend(B.begin(), B.end(), 0, mpl::tag_t{3}));

    requests.waitall();

    int error = -1;
    worker.recv(error, 0);
    if(0 != error) return SUANPAN_FAIL;

    worker.recv(X.begin(), X.end(), 0);

    return SUANPAN_SUCCESS;
}

#endif

#endif

//! @}
