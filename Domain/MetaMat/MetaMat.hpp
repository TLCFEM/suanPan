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
 * @class MetaMat
 * @brief A MetaMat class that holds matrices.
 *
 * @author tlc
 * @date 08/09/2017
 * @version 0.2.0
 * @file MetaMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef METAMAT_HPP
#define METAMAT_HPP

#include "SolverSetting.hpp"
#include "triplet_form.hpp"

#ifdef SUANSPAN_64BIT_INT
using la_it = std::int64_t;
#else
using la_it = std::int32_t;
#endif

template<sp_d T> class MetaMat;

template<sp_d T> class op_add {
    friend MetaMat<T>;

    shared_ptr<MetaMat<T>> mat_a, mat_b;

public:
    explicit op_add(const shared_ptr<MetaMat<T>>& A)
        : mat_a(A)
        , mat_b(nullptr) {}

    op_add(const shared_ptr<MetaMat<T>>& A, const shared_ptr<MetaMat<T>>& B)
        : mat_a(A)
        , mat_b(B) {}
};

template<sp_d T> class op_scale {
    friend MetaMat<T>;

    T scalar;
    op_add<T> bracket;

public:
    op_scale(const T A, const shared_ptr<MetaMat<T>>& B)
        : scalar(A)
        , bracket(B) {}

    op_scale(const T A, op_add<T>&& B)
        : scalar(A)
        , bracket(std::move(B)) {}
};

template<sp_d T> class MetaMat {
protected:
    bool factored = false;

    SolverSetting<T> setting{};

    virtual int direct_solve(Mat<T>&, const Mat<T>&) = 0;

    virtual int direct_solve(Mat<T>&, Mat<T>&&) = 0;

    int direct_solve(Mat<T>& X, const SpMat<T>& B) { return this->direct_solve(X, Mat<T>(B)); }

    int direct_solve(Mat<T>& X, SpMat<T>&& B) { return this->direct_solve(X, B); }

    template<std::invocable<fmat&> F> int mixed_trs(mat& X, mat&& B, F&& trs) {
        auto INFO = 0;

        X = arma::zeros(size(B));

        std::uint8_t counter{0};
        while(counter++ < this->setting.iterative_refinement) {
            const auto multiplier = norm(B);
            if(multiplier < this->setting.tolerance) break;
            suanpan_debug("Mixed precision algorithm multiplier: {:.5E}.\n", multiplier);

            auto residual = conv_to<fmat>::from(B / multiplier);

            if(0 != (INFO = trs(residual))) break;

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;
            B -= this->operator*(incre);
        }

        return INFO;
    }

public:
    triplet_form<T, uword> triplet_mat;

    const uword n_rows;
    const uword n_cols;
    const uword n_elem;

    MetaMat(const uword in_rows, const uword in_cols, const uword in_elem)
        : triplet_mat(in_rows, in_cols)
        , n_rows(in_rows)
        , n_cols(in_cols)
        , n_elem(in_elem) {}

    MetaMat(const MetaMat&) = default;
    MetaMat(MetaMat&&) = delete;
    MetaMat& operator=(const MetaMat&) = delete;
    MetaMat& operator=(MetaMat&&) = delete;
    virtual ~MetaMat() = default;

    void set_solver_setting(const SolverSetting<T>& SS) { setting = SS; }

    [[nodiscard]] SolverSetting<T>& get_solver_setting() { return setting; }

    void set_factored(const bool F) { factored = F; }

    [[nodiscard]] virtual bool is_empty() const = 0;
    virtual void zeros() = 0;

    virtual unique_ptr<MetaMat> make_copy() = 0;

    void unify(const uword K) {
        this->nullify(K);
        this->at(K, K) = T(1);
    }

    virtual void nullify(uword) = 0;

    [[nodiscard]] virtual T max() const = 0;

    /**
     * \brief Access element (read-only), returns zero if out-of-bound
     * \return value
     */
    virtual T operator()(uword, uword) const = 0;
    /**
     * \brief Access element without bound check
     * \return value
     */
    virtual T& unsafe_at(const uword I, const uword J) { return this->at(I, J); }

    /**
     * \brief Access element with bound check
     * \return value
     */
    virtual T& at(uword, uword) = 0;

    [[nodiscard]] virtual const T* memptr() const = 0;
    virtual T* memptr() = 0;

    virtual void scale_accu(T, const shared_ptr<MetaMat>&) = 0;
    virtual void scale_accu(T, const triplet_form<T, uword>&) = 0;

    void operator+=(const shared_ptr<MetaMat>& M) { return this->scale_accu(1., M); }

    void operator-=(const shared_ptr<MetaMat>& M) { return this->scale_accu(-1., M); }

    void operator+=(const op_scale<T>& M) {
        const auto& bracket = M.bracket;
        if(nullptr != bracket.mat_a) this->scale_accu(M.scalar, bracket.mat_a);
        if(nullptr != bracket.mat_b) this->scale_accu(M.scalar, bracket.mat_b);
    }

    void operator-=(const op_scale<T>& M) {
        const auto& bracket = M.bracket;
        if(nullptr != bracket.mat_a) this->scale_accu(-M.scalar, bracket.mat_a);
        if(nullptr != bracket.mat_b) this->scale_accu(-M.scalar, bracket.mat_b);
    }

    void operator+=(const triplet_form<T, uword>& M) { return this->scale_accu(1., M); }

    void operator-=(const triplet_form<T, uword>& M) { return this->scale_accu(-1., M); }

    virtual Mat<T> operator*(const Mat<T>&) const = 0;

    virtual void operator*=(T) = 0;

    template<typename C> requires is_arma_mat<T, C>
    int solve(Mat<T>& X, C&& B) { return this->direct_solve(X, std::forward<C>(B)); }

    template<typename C> requires is_arma_mat<T, C>
    Mat<T> solve(C&& B) {
        Mat<T> X;

        if(SUANPAN_SUCCESS != this->solve(X, std::forward<C>(B))) throw std::runtime_error("fail to solve the system");

        return X;
    }

    [[nodiscard]] virtual int sign_det() const { throw std::runtime_error("not supported"); }

    virtual void allreduce() = 0;

    void save(const char* name) {
        if(!to_mat(*this).save(name, raw_ascii))
            suanpan_error("Cannot save to file \"{}\".\n", name);
    }

    virtual void csc_condense() {}

    virtual void csr_condense() {}
};

template<sp_d T> Mat<T> to_mat(const MetaMat<T>& in_mat) {
    Mat<T> out_mat(in_mat.n_rows, in_mat.n_cols);
    for(uword J = 0; J < in_mat.n_cols; ++J)
        for(uword I = 0; I < in_mat.n_rows; ++I) out_mat(I, J) = in_mat(I, J);
    return out_mat;
}

template<sp_d T> Mat<T> to_mat(const shared_ptr<MetaMat<T>>& in_mat) { return to_mat(*in_mat); }

template<sp_d data_t, sp_i index_t> Mat<data_t> to_mat(const triplet_form<data_t, index_t>& in_mat) {
    Mat<data_t> out_mat(in_mat.n_rows, in_mat.n_cols, fill::zeros);
    for(index_t I = 0; I < in_mat.n_elem; ++I) out_mat(in_mat.row(I), in_mat.col(I)) += in_mat.val(I);
    return out_mat;
}

template<sp_d data_t, sp_i index_t> Mat<data_t> to_mat(const csr_form<data_t, index_t>& in_mat) {
    Mat<data_t> out_mat(in_mat.n_rows, in_mat.n_cols, fill::zeros);

    index_t c_idx = 1;
    for(index_t I = 0; I < in_mat.n_elem; ++I) {
        if(I >= in_mat.row_mem()[c_idx]) ++c_idx;
        out_mat(c_idx - 1, in_mat.col_mem()[I]) += in_mat.val_mem()[I];
    }

    return out_mat;
}

template<sp_d data_t, sp_i index_t> Mat<data_t> to_mat(const csc_form<data_t, index_t>& in_mat) {
    Mat<data_t> out_mat(in_mat.n_rows, in_mat.n_cols, fill::zeros);

    index_t c_idx = 1;
    for(index_t I = 0; I < in_mat.n_elem; ++I) {
        if(I >= in_mat.col_mem()[c_idx]) ++c_idx;
        out_mat(in_mat.row_mem()[I], c_idx - 1) += in_mat.val_mem()[I];
    }

    return out_mat;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> to_triplet_form(MetaMat<data_t>* in_mat) {
    if(!in_mat->triplet_mat.is_empty()) return triplet_form<data_t, index_t>(in_mat->triplet_mat);

    const sp_i auto n_rows = index_t(in_mat->n_rows);
    const sp_i auto n_cols = index_t(in_mat->n_cols);
    const sp_i auto n_elem = index_t(in_mat->n_elem);

    triplet_form<data_t, index_t> out_mat(n_rows, n_cols, n_elem);
    for(index_t J = 0; J < n_cols; ++J)
        for(index_t I = 0; I < n_rows; ++I) out_mat.at(I, J) = in_mat->operator()(I, J);

    return out_mat;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> to_triplet_form(const shared_ptr<MetaMat<data_t>>& in_mat) { return to_triplet_form<data_t, index_t>(in_mat.get()); }

#endif

//! @}
