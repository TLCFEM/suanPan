/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class SparseMatMPI
 * @brief A SparseMatMPI class that holds matrices.
 *
 * * MUMPS uses int.
 *
 * @author tlc
 * @date 14/06/2023
 * @version 0.1.0
 * @file SparseMatMPI.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppClangTidyClangDiagnosticMissingFieldInitializers
#ifndef SPARSEMATMPI_HPP
#define SPARSEMATMPI_HPP

#if defined(SUANPAN_MPI) && defined(SUANPAN_MKL)

#include <mpi.h>

extern int SUANPAN_NUM_NODES;

template<sp_d T> class SparseMatMPIPARDISO final : public SparseMat<T> {
    int iparm[64];

protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatMPIPARDISO(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem)
        , iparm{} {
        iparm[34] = 1; // zero-based indexing
    }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatMPIPARDISO>(*this); }
};

template<sp_d T> int SparseMatMPIPARDISO<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    X.set_size(B.n_rows, B.n_cols);

    csr_form<T, int> csr_mat(this->triplet_mat, SparseBase::ZERO, true);

    const auto n = static_cast<int>(B.n_rows);
    const auto nrhs = static_cast<int>(B.n_cols);
    const auto nnz = static_cast<int>(csr_mat.n_elem);

    MPI_Comm worker;
    MPI_Comm_spawn("solver.pardiso", MPI_ARGV_NULL, SUANPAN_NUM_NODES, MPI_INFO_NULL, 0, MPI_COMM_SELF, &worker, MPI_ERRCODES_IGNORE);

    int config[8];

    config[0] = 11;   // mtype
    config[1] = nrhs; // nrhs
    config[2] = 1;    // maxfct
    config[3] = 1;    // mnum
    config[4] = 0;    // msglvl
    config[5] = n;    // n
    config[6] = nnz;  // nnz
    config[7] = std::is_same_v<T, double> ? 1 : -1;
    const auto FLOAT_TYPE = std::is_same_v<T, double> ? MPI_DOUBLE : MPI_FLOAT;

    MPI_Comm remote;
    MPI_Intercomm_merge(worker, 0, &remote);
    MPI_Bcast(&config, 8, MPI_INT, 0, remote);

    MPI_Request requests[5];
    MPI_Isend(&iparm, 64, MPI_INT, 0, 0, worker, &requests[0]);
    MPI_Isend(csr_mat.row_mem(), n + 1, MPI_INT, 0, 0, worker, &requests[1]);
    MPI_Isend(csr_mat.col_mem(), nnz, MPI_INT, 0, 0, worker, &requests[2]);
    MPI_Isend(csr_mat.val_mem(), nnz, FLOAT_TYPE, 0, 0, worker, &requests[3]);
    MPI_Isend(B.memptr(), static_cast<int>(B.n_elem), FLOAT_TYPE, 0, 0, worker, &requests[4]);
    MPI_Waitall(5, requests, MPI_STATUSES_IGNORE);

    int error = -1;
    MPI_Recv(&error, 1, MPI_INT, 0, 0, worker, MPI_STATUS_IGNORE);
    if(0 == error) {
        MPI_Recv(X.memptr(), static_cast<int>(B.n_elem), FLOAT_TYPE, 0, 0, worker, MPI_STATUS_IGNORE);
        return SUANPAN_SUCCESS;
    }

    return SUANPAN_FAIL;
}

#endif

#endif

//! @}
