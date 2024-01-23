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

#ifndef OPERATOR_TIMES_HPP
#define OPERATOR_TIMES_HPP

template<sp_d T> op_add<T> operator+(const shared_ptr<MetaMat<T>>& A, const shared_ptr<MetaMat<T>>& B) { return op_add<T>(A, B); }

template<sp_d T> op_scale<T> operator*(const T value, const shared_ptr<MetaMat<T>>& M) { return op_scale<T>(value, M); }

template<sp_d T> op_scale<T> operator*(const T value, op_add<T>&& M) { return op_scale<T>(value, std::forward<op_add<T>>(M)); }

template<sp_d T> const shared_ptr<MetaMat<T>>& operator+=(const shared_ptr<MetaMat<T>>& M, const op_scale<T>& A) {
    M->operator+=(A);
    return M;
}

template<sp_d T> const unique_ptr<MetaMat<T>>& operator+=(const unique_ptr<MetaMat<T>>& M, const op_scale<T>& A) {
    M->operator+=(A);
    return M;
}

template<sp_d T> unique_ptr<MetaMat<T>> operator*(const T value, unique_ptr<MetaMat<T>>&& M) {
    if(nullptr == M) return nullptr;

    M->operator*=(value);
    return std::forward<unique_ptr<MetaMat<T>>>(M);
}

//template<sp_d T> unique_ptr<MetaMat<T>> operator*(const T value, const shared_ptr<MetaMat<T>>& M) {
//    if(nullptr == M) return nullptr;
//
//    auto N = M->make_copy();
//    N->operator*=(value);
//    return N;
//}

template<sp_d T> Mat<T> operator*(const shared_ptr<MetaMat<T>>& M, const Mat<T>& A) { return M->operator*(A); }

template<sp_d T> Mat<T> operator*(const unique_ptr<MetaMat<T>>& M, const Mat<T>& A) { return M->operator*(A); }

template<sp_d T> const shared_ptr<MetaMat<T>>& operator*=(const shared_ptr<MetaMat<T>>& M, const T value) {
    M->operator*=(value);
    return M;
}

//template<sp_d T> unique_ptr<MetaMat<T>> operator+(const shared_ptr<MetaMat<T>>& A, const shared_ptr<MetaMat<T>>& B) {
//    auto M = B->make_copy();
//    M->operator+=(A);
//    return M;
//}

template<sp_d T> unique_ptr<MetaMat<T>> operator+(const shared_ptr<MetaMat<T>>& A, unique_ptr<MetaMat<T>>&& B) {
    B->operator+=(A);
    return std::forward<unique_ptr<MetaMat<T>>>(B);
}

template<sp_d T> unique_ptr<MetaMat<T>> operator+(unique_ptr<MetaMat<T>>&& A, unique_ptr<MetaMat<T>>&& B) {
    A->operator+=(std::forward<unique_ptr<MetaMat<T>>>(B));
    return std::forward<unique_ptr<MetaMat<T>>>(A);
}

template<sp_d T> const shared_ptr<MetaMat<T>>& operator+=(const shared_ptr<MetaMat<T>>& M, const shared_ptr<MetaMat<T>>& A) {
    M->operator+=(A);
    return M;
}

template<sp_d T> const shared_ptr<MetaMat<T>>& operator+=(const shared_ptr<MetaMat<T>>& M, unique_ptr<MetaMat<T>>&& A) {
    M->operator+=(std::forward<unique_ptr<MetaMat<T>>>(A));
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

#endif // OPERATOR_TIMES_HPP
