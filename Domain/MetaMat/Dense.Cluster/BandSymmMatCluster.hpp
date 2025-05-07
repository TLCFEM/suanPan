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
 * @class BandSymmMatCluster
 * @brief A BandSymmMatCluster class that holds matrices.
 *
 * @author tlc
 * @date 31/03/2025
 * @version 0.1.0
 * @file BandSymmMatCluster.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDSYMMMATCLUSTER_HPP
#define BANDSYMMMATCLUSTER_HPP

#include "../DenseMat.hpp"

#include <ezp/ezp/ppbsv.hpp>

template<sp_d T> class BandSymmMatCluster final : public DenseMat<T> {
    static constexpr char UPLO = 'L';

    using solver_t = ezp::ppbsv<T, la_it, UPLO>;
    using indexer_t = typename solver_t::indexer;

    static T bin;

    solver_t solver;
    indexer_t indexer;

    const uword band;

    int solve_trs(Mat<T>&, Mat<T>&&);

protected:
    using DenseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    BandSymmMatCluster(const uword in_size, const uword in_bandwidth)
        : DenseMat<T>(in_size, in_size, (in_bandwidth + 1) * in_size)
        , solver()
        , indexer(in_size, in_bandwidth)
        , band(in_bandwidth) {}

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<BandSymmMatCluster>(*this); }

    void nullify(const uword K) override {
        this->factored = false;
        suanpan::for_each(std::max(K, band) - band, K, [&](const uword I) { this->memory[indexer(K, I)] = T(0); });
        suanpan::for_each(K, std::min(this->n_rows, K + band + 1), [&](const uword I) { this->memory[indexer(I, K)] = T(0); });
    }

    [[nodiscard]] SpMat<T> extract_col(const uword K) override {
        SpMat<T> output(this->n_rows, 1);
        for(auto I = std::max(K, band) - band; I < K; ++I) output.at(I, 0) = this->memory[indexer(K, I)];
        for(auto I = K; I < std::min(this->n_rows, K + band + 1); ++I) output.at(I, 0) = this->memory[indexer(I, K)];
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

template<sp_d T> T BandSymmMatCluster<T>::bin = T(0);

template<sp_d T> Mat<T> BandSymmMatCluster<T>::operator*(const Mat<T>& X) const {
    static constexpr blas_int INC = 1;
    static constexpr T ALPHA = T(1), BETA = T(0);

    Mat<T> Y(arma::size(X));

    const auto N = static_cast<blas_int>(this->n_cols);
    const auto K = static_cast<blas_int>(band);
    const auto LDA = K + 1;

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_ssbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }
    else {
        using E = double;
        suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_dsbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }

    return Y;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnarrowing"
template<sp_d T> int BandSymmMatCluster<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::move(B));

    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw std::invalid_argument("requires a square matrix"); });

    this->factored = true;

    const auto INFO = bcast_from_root(solver.solve({this->n_rows, this->n_cols, this->band, this->memptr()}, {B.n_rows, B.n_cols, B.memptr()}));

    if(0 == INFO) bcast_from_root(X = std::move(B));
    else suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int BandSymmMatCluster<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    const auto INFO = bcast_from_root(solver.solve({B.n_rows, B.n_cols, B.memptr()}));

    if(0 == INFO) bcast_from_root(X = std::move(B));
    else suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}
#pragma GCC diagnostic pop

#endif

//! @}
