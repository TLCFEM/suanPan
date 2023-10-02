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

#ifndef OPERATOR_TIMES_HPP
#define OPERATOR_TIMES_HPP

#include "FullMat.hpp"

template<sp_d T> unique_ptr<MetaMat<T>> operator*(const T value, const unique_ptr<MetaMat<T>>& M) {
    if(nullptr == M) return nullptr;

    auto N = M->make_copy();
    N->operator*=(value);
    return N;
}

template<sp_d T> unique_ptr<MetaMat<T>> operator*(const T value, const shared_ptr<MetaMat<T>>& M) {
    if(nullptr == M) return nullptr;

    auto N = M->make_copy();
    N->operator*=(value);
    return N;
}

template<sp_d T> const shared_ptr<MetaMat<T>>& operator*=(const shared_ptr<MetaMat<T>>& M, const T value) {
    M->operator*=(value);
    return M;
}

template<sp_d T> const unique_ptr<MetaMat<T>>& operator*=(const unique_ptr<MetaMat<T>>& M, const T value) {
    M->operator*=(value);
    return M;
}

template<sp_d T> const shared_ptr<MetaMat<T>>& operator+=(const shared_ptr<MetaMat<T>>& M, const shared_ptr<MetaMat<T>>& A) {
    M->operator+=(A);
    return M;
}

template<sp_d T> const unique_ptr<MetaMat<T>>& operator+=(const unique_ptr<MetaMat<T>>& M, const shared_ptr<MetaMat<T>>& A) {
    M->operator+=(A);
    return M;
}

template<sp_d DT, sp_i IT> const unique_ptr<MetaMat<DT>>& operator+=(const unique_ptr<MetaMat<DT>>& M, const triplet_form<DT, IT>& A) {
    M->operator+=(A);
    return M;
}

template<sp_d T> const shared_ptr<MetaMat<T>>& operator+=(const shared_ptr<MetaMat<T>>& M, unique_ptr<MetaMat<T>>&& A) {
    M->operator+=(std::forward<unique_ptr<MetaMat<T>>>(A));
    return M;
}

template<sp_d T> unique_ptr<MetaMat<T>> operator+(unique_ptr<MetaMat<T>>&& A, unique_ptr<MetaMat<T>>&& B) {
    if(nullptr == A && nullptr == B) return nullptr;

    if(nullptr != A) {
        A->operator+=(std::forward<unique_ptr<MetaMat<T>>>(B));
        return std::forward<unique_ptr<MetaMat<T>>>(A);
    }

    return std::forward<unique_ptr<MetaMat<T>>>(B);
}

template<sp_d T> unique_ptr<MetaMat<T>> operator+(const shared_ptr<MetaMat<T>>& A, unique_ptr<MetaMat<T>>&& B) {
    B->operator+=(A);
    return std::forward<unique_ptr<MetaMat<T>>>(B);
}

template<sp_d T> unique_ptr<MetaMat<T>> operator+(unique_ptr<MetaMat<T>>&& B, const shared_ptr<MetaMat<T>>& A) {
    B->operator+=(A);
    return std::forward<unique_ptr<MetaMat<T>>>(B);
}

template<sp_d T> const shared_ptr<MetaMat<T>>& operator+(const shared_ptr<MetaMat<T>>& A, const shared_ptr<MetaMat<T>>& B) {
    B->operator+=(A);
    return B;
}

template<sp_d T> const shared_ptr<MetaMat<T>>& operator-=(const shared_ptr<MetaMat<T>>& M, const shared_ptr<MetaMat<T>>& A) {
    M->operator-=(A);
    return M;
}

template<sp_d T> const shared_ptr<MetaMat<T>>& operator-=(const shared_ptr<MetaMat<T>>& M, unique_ptr<MetaMat<T>>&& A) {
    M->operator-=(std::forward<unique_ptr<MetaMat<T>>>(A));
    return M;
}

template<sp_d T> unique_ptr<MetaMat<T>> operator-(unique_ptr<MetaMat<T>>&& A, unique_ptr<MetaMat<T>>&& B) {
    if(nullptr == A && nullptr == B) return nullptr;

    if(nullptr != A) {
        A->operator-=(std::forward<unique_ptr<MetaMat<T>>>(B));
        return std::forward<unique_ptr<MetaMat<T>>>(A);
    }

    return std::forward<unique_ptr<MetaMat<T>>>(B);
}

template<sp_d T> unique_ptr<MetaMat<T>> operator-(const shared_ptr<MetaMat<T>>& A, unique_ptr<MetaMat<T>>&& B) {
    B->operator-=(A);
    return std::forward<unique_ptr<MetaMat<T>>>(B);
}

template<sp_d T> unique_ptr<MetaMat<T>> operator-(unique_ptr<MetaMat<T>>&& B, const shared_ptr<MetaMat<T>>& A) {
    B->operator-=(A);
    return std::forward<unique_ptr<MetaMat<T>>>(B);
}

template<sp_d T> Mat<T> operator*(const shared_ptr<MetaMat<T>>& M, const Mat<T>& A) {
    if(nullptr == M) return nullptr;

    return M->operator*(A);
}

template<sp_d T> Mat<T> operator*(const unique_ptr<MetaMat<T>>& M, const Mat<T>& A) {
    if(nullptr == M) return nullptr;

    return M->operator*(A);
}

template<sp_d T> Mat<T> operator*(const Mat<T>& A, const FullMat<T>& B) {
    Mat<T> C(A.n_rows, A.n_cols);

    constexpr auto TRAN = 'N';

    const auto M = static_cast<int>(A.n_rows);
    const auto N = static_cast<int>(B.n_cols);
    const auto K = static_cast<int>(A.n_cols);
    T ALPHA = 1.;
    const auto LDA = M;
    const auto LDB = K;
    T BETA = 0.;
    const auto LDC = M;

    if(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_sgemm)(&TRAN, &TRAN, &M, &N, &K, (E*)&ALPHA, (E*)A.memptr(), &LDA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
    }
    else if(std::is_same_v<T, double>) {
        using E = double;
        arma_fortran(arma_dgemm)(&TRAN, &TRAN, &M, &N, &K, (E*)&ALPHA, (E*)A.memptr(), &LDA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
    }

    return C;
}

template<sp_d T, sp_i IT> triplet_form<T, IT> operator*(const T value, const triplet_form<T, IT>& M) {
    auto N = M;
    N *= value;
    return N;
}

template<sp_d T, sp_i IT> triplet_form<T, IT> operator*(const T value, triplet_form<T, IT>&& M) {
    M *= value;
    return M;
}

#endif // OPERATOR_TIMES_HPP
