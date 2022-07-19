/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "triplet_form.hpp"
#include "IterativeSolver.hpp"

template<sp_d T> class MetaMat {
protected:
    bool factored = false;

    SolverSetting<T> setting{};

public:
    triplet_form<T, uword> triplet_mat;

    const uword n_rows;
    const uword n_cols;
    const uword n_elem;

    MetaMat(uword, uword, uword);
    MetaMat(const MetaMat&) = default;
    MetaMat(MetaMat&&) noexcept = delete;
    MetaMat& operator=(const MetaMat&) = delete;
    MetaMat& operator=(MetaMat&&) noexcept = delete;
    virtual ~MetaMat() = default;

    void set_solver_setting(const SolverSetting<T>&);
    [[nodiscard]] SolverSetting<T>& get_solver_setting();

    void set_factored(bool);

    [[nodiscard]] virtual bool is_empty() const = 0;
    virtual void zeros() = 0;

    virtual unique_ptr<MetaMat> make_copy() = 0;

    virtual void unify(uword) = 0;
    virtual void nullify(uword) = 0;

    [[nodiscard]] virtual T max() const = 0;
    [[nodiscard]] virtual Col<T> diag() const = 0;

    virtual const T& operator()(uword, uword) const = 0;
    virtual T& at(uword, uword) = 0;

    [[nodiscard]] virtual const T* memptr() const = 0;
    virtual T* memptr() = 0;

    virtual void operator+=(const shared_ptr<MetaMat>&) = 0;
    virtual void operator-=(const shared_ptr<MetaMat>&) = 0;

    virtual void operator+=(const triplet_form<T, uword>&) = 0;
    virtual void operator-=(const triplet_form<T, uword>&) = 0;

    virtual Mat<T> operator*(const Mat<T>&) const = 0;

    virtual void operator*=(T) = 0;

    Mat<T> solve(const Mat<T>&);
    Mat<T> solve(const SpMat<T>&);
    Mat<T> solve(Mat<T>&&);
    Mat<T> solve(SpMat<T>&&);

    virtual int solve(Mat<T>&, const Mat<T>&) = 0;
    virtual int solve(Mat<T>&, const SpMat<T>&);
    virtual int solve(Mat<T>&, Mat<T>&&);
    virtual int solve(Mat<T>&, SpMat<T>&&);

    [[nodiscard]] virtual int sign_det() const = 0;

    void save(const char*);

    virtual void csc_condense();
    virtual void csr_condense();

    int iterative_solve(Mat<T>&, const Mat<T>&);

    [[nodiscard]] Col<T> evaluate(const Col<T>&) const;
};

template<sp_d T> MetaMat<T>::MetaMat(const uword in_rows, const uword in_cols, const uword in_elem)
    : triplet_mat(in_rows, in_cols)
    , n_rows(in_rows)
    , n_cols(in_cols)
    , n_elem(in_elem) {}

template<sp_d T> void MetaMat<T>::set_solver_setting(const SolverSetting<T>& SS) { setting = SS; }

template<sp_d T> SolverSetting<T>& MetaMat<T>::get_solver_setting() { return setting; }

template<sp_d T> void MetaMat<T>::set_factored(const bool F) { factored = F; }

template<sp_d T> Mat<T> MetaMat<T>::solve(const Mat<T>& B) {
    Mat<T> X;
    if(0 != this->solve(X, B)) X.reset();
    return X;
}

template<sp_d T> Mat<T> MetaMat<T>::solve(const SpMat<T>& B) {
    Mat<T> X;
    if(0 != this->solve(X, B)) X.reset();
    return X;
}

template<sp_d T> Mat<T> MetaMat<T>::solve(Mat<T>&& B) {
    Mat<T> X;
    if(0 != this->solve(X, std::forward<Mat<T>>(B))) X.reset();
    return X;
}

template<sp_d T> Mat<T> MetaMat<T>::solve(SpMat<T>&& B) {
    Mat<T> X;
    if(0 != this->solve(X, std::forward<SpMat<T>>(B))) X.reset();
    return X;
}

template<sp_d T> int MetaMat<T>::solve(Mat<T>& X, const SpMat<T>& B) { return this->solve(X, Mat<T>(B)); }

template<sp_d T> int MetaMat<T>::solve(Mat<T>& X, Mat<T>&& B) { return this->solve(X, B); }

template<sp_d T> int MetaMat<T>::solve(Mat<T>& X, SpMat<T>&& B) { return this->solve(X, B); }

template<sp_d T> void MetaMat<T>::save(const char* name) { if(!to_mat(*this).save(name)) suanpan_error("cannot save matrix to file.\n"); }

template<sp_d T> void MetaMat<T>::csc_condense() {}

template<sp_d T> void MetaMat<T>::csr_condense() {}

template<sp_d T> int MetaMat<T>::iterative_solve(Mat<T>& X, const Mat<T>& B) {
    X = B;

    SimpleJacobi<T> preconditioner(this->diag());

    std::atomic_int code = 0;

    if(IterativeSolver::GMRES == setting.iterative_solver)
        suanpan_for(0llu, B.n_cols, [&](const uword I) {
            Col<T> sub_x(X.colptr(I), X.n_rows, false, true);
            const Col<T> sub_b(B.colptr(I), B.n_rows);
            auto col_setting = setting;
            code += GMRES(this, sub_x, sub_b, preconditioner, col_setting);
        });
    else if(IterativeSolver::BICGSTAB == setting.iterative_solver)
        suanpan_for(0llu, B.n_cols, [&](const uword I) {
            Col<T> sub_x(X.colptr(I), X.n_rows, false, true);
            const Col<T> sub_b(B.colptr(I), B.n_rows);
            auto col_setting = setting;
            code += BiCGSTAB(this, sub_x, sub_b, preconditioner, col_setting);
        });
    else throw invalid_argument("no proper iterative solver assigned but somehow iterative solving is called");

    return 0 == code ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

template<sp_d T> Col<T> MetaMat<T>::evaluate(const Col<T>& X) const { return this->operator*(X); }

template<sp_d T> Mat<T> to_mat(const MetaMat<T>& in_mat) {
    Mat<T> out_mat(in_mat.n_rows, in_mat.n_cols);
    for(uword J = 0; J < in_mat.n_cols; ++J) for(uword I = 0; I < in_mat.n_rows; ++I) out_mat(I, J) = in_mat(I, J);
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

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> to_triplet_form(const shared_ptr<MetaMat<data_t>>& in_mat) {
    if(!in_mat->triplet_mat.is_empty()) return triplet_form<data_t, index_t>(in_mat->triplet_mat);

    const sp_i auto n_rows = index_t(in_mat->n_rows);
    const sp_i auto n_cols = index_t(in_mat->n_cols);
    const sp_i auto n_elem = index_t(in_mat->n_elem);

    triplet_form<data_t, index_t> out_mat(n_rows, n_cols, n_elem);
    for(index_t J = 0; J < n_cols; ++J) for(index_t I = 0; I < n_rows; ++I) out_mat.at(I, J) = in_mat->operator()(I, J);

    return out_mat;
}

#endif

//! @}
