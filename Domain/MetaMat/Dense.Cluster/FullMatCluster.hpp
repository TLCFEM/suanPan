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
 * @class FullMatCluster
 * @brief A FullMatCluster class that holds matrices.
 *
 * @author tlc
 * @date 31/03/2025
 * @version 0.1.0
 * @file FullMatCluster.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef FULLMATCLUSTER_HPP
#define FULLMATCLUSTER_HPP

#include "../DenseMat.hpp"

#include <ezp/ezp/pgesv.hpp>
#include <ezp/ezp/pposv.hpp>

template<sp_d T, typename solver_t> class FullMatBaseCluster : public DenseMat<T> {
    solver_t solver;

    int solve_trs(Mat<T>&, Mat<T>&&);

protected:
    using DenseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    FullMatBaseCluster(const uword in_rows, const uword in_cols)
        : DenseMat<T>(in_rows, in_cols, in_rows * in_cols)
        , solver() {}

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<FullMatBaseCluster>(*this); }

    void nullify(const uword K) override {
        this->factored = false;
        suanpan::for_each(this->n_rows, [&](const uword I) { this->at(I, K) = T(0); });
        suanpan::for_each(this->n_cols, [&](const uword I) { this->at(K, I) = T(0); });
    }

    T operator()(const uword in_row, const uword in_col) const override { return this->memory[in_row + in_col * this->n_rows]; }

    T& at(const uword in_row, const uword in_col) override {
        this->factored = false;
        return this->memory[in_row + in_col * this->n_rows];
    }

    Mat<T> operator*(const Mat<T>&) const override;
};

template<sp_d T, typename solver_t> Mat<T> FullMatBaseCluster<T, solver_t>::operator*(const Mat<T>& B) const {
    static constexpr char TRAN = 'N';
    static constexpr T ALPHA = T(1), BETA = T(0);

    Mat<T> C(arma::size(B));

    const auto M = static_cast<blas_int>(this->n_rows);
    const auto N = static_cast<blas_int>(this->n_cols);

    if(1 == B.n_cols) {
        static constexpr blas_int INC = 1;

        if constexpr(std::is_same_v<T, float>) {
            using E = float;
            arma_fortran(arma_sgemv)(&TRAN, &M, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &INC, (E*)&BETA, (E*)C.memptr(), &INC);
        }
        else {
            using E = double;
            arma_fortran(arma_dgemv)(&TRAN, &M, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &INC, (E*)&BETA, (E*)C.memptr(), &INC);
        }
    }
    else {
        const auto K = static_cast<blas_int>(B.n_cols);

        if constexpr(std::is_same_v<T, float>) {
            using E = float;
            arma_fortran(arma_sgemm)(&TRAN, &TRAN, &M, &K, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &N, (E*)&BETA, (E*)C.memptr(), &M);
        }
        else {
            using E = double;
            arma_fortran(arma_dgemm)(&TRAN, &TRAN, &M, &K, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &N, (E*)&BETA, (E*)C.memptr(), &M);
        }
    }

    return C;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnarrowing"
template<sp_d T, typename solver_t> int FullMatBaseCluster<T, solver_t>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::move(B));

    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw std::invalid_argument("requires a square matrix"); });

    this->factored = true;

    const auto INFO = bcast_from_root(solver.solve({this->n_rows, this->n_cols, this->memptr()}, {B.n_rows, B.n_cols, B.memptr()}));

    if(0 == INFO) bcast_from_root(X = std::move(B));
    else suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T, typename solver_t> int FullMatBaseCluster<T, solver_t>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    const auto INFO = bcast_from_root(solver.solve({B.n_rows, B.n_cols, B.memptr()}));

    if(0 == INFO) bcast_from_root(X = std::move(B));
    else suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}
#pragma GCC diagnostic pop

template<sp_d T> using FullMatCluster = FullMatBaseCluster<T, ezp::pgesv<T, la_it>>;
template<sp_d T> using FullSymmMatCluster = FullMatBaseCluster<T, ezp::pposv<T, la_it>>;

#endif

//! @}
