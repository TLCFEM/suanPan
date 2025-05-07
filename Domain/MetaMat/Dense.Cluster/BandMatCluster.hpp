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
 * @class BandMatCluster
 * @brief A BandMatCluster class that holds matrices.
 *
 * @author tlc
 * @date 31/03/2025
 * @version 0.1.0
 * @file BandMatCluster.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDMATCLUSTER_HPP
#define BANDMATCLUSTER_HPP

#include "../DenseMat.hpp"

#include <ezp/ezp/pgbsv.hpp>

template<sp_d T> class BandMatCluster : public DenseMat<T> {
    using solver_t = ezp::pgbsv<T, la_it>;
    using indexer_t = typename solver_t::indexer;

    static T bin;

    solver_t solver;
    indexer_t indexer;

    int solve_trs(Mat<T>&, Mat<T>&&);

protected:
    const uword l_band, u_band;

    using DenseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    BandMatCluster(const uword in_size, const uword in_l, const uword in_u)
        : DenseMat<T>(in_size, in_size, (2 * (in_l + in_u) + 1) * in_size)
        , solver()
        , indexer(in_size, in_l, in_u)
        , l_band(in_l)
        , u_band(in_u) {
        if(2 * (in_l + in_u) + 1 >= in_size)
            suanpan_warning("The storage requirement for the banded matrix is larger than that of a full matrix, consider using a full/sparse matrix instead.\n");
    }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<BandMatCluster>(*this); }

    void nullify(const uword K) override {
        this->factored = false;
        suanpan::for_each(std::max(K, u_band) - u_band, std::min(this->n_rows, K + l_band + 1), [&](const uword I) { this->memory[2 * u_band + l_band + I + 2 * K * (l_band + u_band)] = T(0); });
        suanpan::for_each(std::max(K, l_band) - l_band, std::min(this->n_cols, K + u_band + 1), [&](const uword I) { this->memory[2 * u_band + l_band + K + 2 * I * (l_band + u_band)] = T(0); });
    }

    [[nodiscard]] SpMat<T> extract_col(const uword K) override {
        SpMat<T> output(this->n_rows, 1);
        for(auto I = std::max(K, u_band) - u_band; I < std::min(this->n_rows, K + l_band + 1); ++I) output.at(I, 0) = this->memory[2 * u_band + l_band + I + 2 * K * (l_band + u_band)];
        return output;
    }

    T operator()(const uword in_row, const uword in_col) const override {
        const auto pos = indexer(in_row, in_col);
        if(pos < 0) [[unlikely]]
            return bin = T(0);
        return this->memory[pos];
    }

    T& unsafe_at(const uword in_row, const uword in_col) override {
        this->factored = false;
        return this->memory[indexer(in_row, in_col)];
    }

    T& at(const uword in_row, const uword in_col) override {
        const auto pos = indexer(in_row, in_col);
        if(pos < 0) [[unlikely]]
            return bin = T(0);
        this->factored = false;
        return this->memory[pos];
    }

    Mat<T> operator*(const Mat<T>&) const override;
};

template<sp_d T> T BandMatCluster<T>::bin = T(0);

template<sp_d T> Mat<T> BandMatCluster<T>::operator*(const Mat<T>& X) const {
    static constexpr char TRAN = 'N';
    static constexpr blas_int INC = 1;
    static constexpr T ALPHA = T(1), BETA = T(0);

    Mat<T> Y(arma::size(X));

    const auto s_band = l_band + u_band;

    const auto M = static_cast<blas_int>(this->n_rows);
    const auto N = static_cast<blas_int>(this->n_cols);
    const auto KL = static_cast<blas_int>(l_band);
    const auto KU = static_cast<blas_int>(u_band);
    const auto LDA = static_cast<blas_int>(2 * s_band + 1);

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + s_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }
    else {
        using E = double;
        suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + s_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }

    return Y;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnarrowing"
template<sp_d T> int BandMatCluster<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::move(B));

    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw std::invalid_argument("requires a square matrix"); });

    this->factored = true;

    const auto INFO = bcast_from_root(solver.solve({this->n_rows, this->n_cols, this->l_band, this->u_band, this->memptr()}, {B.n_rows, B.n_cols, B.memptr()}));

    if(0 == INFO) bcast_from_root(X = std::move(B));
    else suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int BandMatCluster<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    const auto INFO = bcast_from_root(solver.solve({B.n_rows, B.n_cols, B.memptr()}));

    if(0 == INFO) bcast_from_root(X = std::move(B));
    else suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}
#pragma GCC diagnostic pop

#endif

//! @}
