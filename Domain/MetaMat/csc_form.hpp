﻿/*******************************************************************************
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

#ifndef CSC_FORM_HPP
#define CSC_FORM_HPP

#include "triplet_form.hpp"

template<sp_d data_t, sp_i index_t> class csr_form;

template<sp_d data_t, sp_i index_t> class csc_form final {
    const data_t bin = data_t(0);

    using index_ptr = std::unique_ptr<index_t[]>;
    using data_ptr = std::unique_ptr<data_t[]>;

    index_ptr row_idx = nullptr; // index storage
    index_ptr col_ptr = nullptr; // index storage
    data_ptr val_idx = nullptr;  // value storage

    template<sp_d in_dt, sp_i in_it> void copy_to(in_it* const new_row_idx, in_it* const new_col_ptr, in_dt* const new_val_idx) const {
        suanpan::for_each(n_cols + 1, [&](const index_t I) { new_col_ptr[I] = in_it(col_ptr[I]); });
        suanpan::for_each(n_elem, [&](const index_t I) {
            new_row_idx[I] = in_it(row_idx[I]);
            new_val_idx[I] = in_dt(val_idx[I]);
        });
    }

    void init(const index_t in_elem) {
        row_idx = std::move(index_ptr(new index_t[in_elem]));
        col_ptr = std::move(index_ptr(new index_t[n_cols + 1]));
        val_idx = std::move(data_ptr(new data_t[in_elem]));
    }

public:
    const index_t n_rows = 0;
    const index_t n_cols = 0;
    const index_t n_elem = 0;

    csc_form() = default;
    csc_form(const csc_form&);
    csc_form(csc_form&&) noexcept;
    csc_form& operator=(const csc_form&);
    csc_form& operator=(csc_form&&) noexcept;
    ~csc_form() = default;

    [[nodiscard]] const index_t* row_mem() const { return row_idx.get(); }

    [[nodiscard]] const index_t* col_mem() const { return col_ptr.get(); }

    [[nodiscard]] const data_t* val_mem() const { return val_idx.get(); }

    [[nodiscard]] index_t* row_mem() { return row_idx.get(); }

    [[nodiscard]] index_t* col_mem() { return col_ptr.get(); }

    [[nodiscard]] data_t* val_mem() { return val_idx.get(); }

    index_t row(const index_t I) const { return row_idx[I]; }

    index_t col(const index_t I) const { return col_ptr[I]; }

    data_t val(const index_t I) const { return val_idx[I]; }

    [[nodiscard]] data_t max() const {
        if(0 == n_elem) return data_t(0);
        return *suanpan::max_element(val_idx.get(), val_idx.get() + n_elem);
    }

    void print() const;

    template<sp_d T2> csc_form operator*(const T2 scalar) const {
        auto copy = *this;
        return copy *= scalar;
    }

    template<sp_d T2> csc_form operator/(const T2 scalar) const {
        auto copy = *this;
        return copy /= scalar;
    }

    template<sp_d T2> csc_form& operator*=(const T2 scalar) {
        suanpan_for_each(val_idx.get(), val_idx.get() + n_elem, [=](data_t& I) { I *= data_t(scalar); });
        return *this;
    }

    template<sp_d T2> csc_form& operator/=(const T2 scalar) {
        suanpan_for_each(val_idx.get(), val_idx.get() + n_elem, [=](data_t& I) { I /= data_t(scalar); });
        return *this;
    }

    template<sp_d in_dt, sp_i in_it> explicit csc_form(triplet_form<in_dt, in_it>&, SparseBase = SparseBase::ZERO, bool = false);

    template<sp_d in_dt, sp_i in_it> csc_form& operator=(triplet_form<in_dt, in_it>&);

    data_t operator()(const index_t in_row, const index_t in_col) const {
        if(in_row < n_rows && in_col < n_cols)
            for(auto I = col_ptr[in_col]; I < col_ptr[in_col + 1]; ++I)
                if(in_row == row_idx[I]) return val_idx[I];
        return access::rw(bin) = data_t(0);
    }

    Mat<data_t> operator*(const Col<data_t>& in_mat) const {
        Mat<data_t> out_mat = arma::zeros<Mat<data_t>>(in_mat.n_rows, 1);
        for(index_t I = 0; I < n_cols; ++I)
            for(auto J = col_ptr[I]; J < col_ptr[I + 1]; ++J) out_mat(row_idx[J]) += val_idx[J] * in_mat(I);
        return out_mat;
    }

    Mat<data_t> operator*(const Mat<data_t>& in_mat) const {
        Mat<data_t> out_mat = arma::zeros<Mat<data_t>>(in_mat.n_rows, in_mat.n_cols);
        for(index_t I = 0; I < n_cols; ++I)
            for(auto J = col_ptr[I]; J < col_ptr[I + 1]; ++J) out_mat.row(row_idx[J]) += val_idx[J] * in_mat.row(I);
        return out_mat;
    }
};

template<sp_d data_t, sp_i index_t> csc_form<data_t, index_t>::csc_form(const csc_form& in_mat)
    : n_rows{in_mat.n_rows}
    , n_cols{in_mat.n_cols}
    , n_elem{in_mat.n_elem} {
    init(n_elem);
    in_mat.copy_to(row_idx.get(), col_ptr.get(), val_idx.get());
}

template<sp_d data_t, sp_i index_t> csc_form<data_t, index_t>::csc_form(csc_form&& in_mat) noexcept
    : row_idx{std::move(in_mat.row_idx)}
    , col_ptr{std::move(in_mat.col_ptr)}
    , val_idx{std::move(in_mat.val_idx)}
    , n_rows{in_mat.n_rows}
    , n_cols{in_mat.n_cols}
    , n_elem{in_mat.n_elem} {}

template<sp_d data_t, sp_i index_t> csc_form<data_t, index_t>& csc_form<data_t, index_t>::operator=(const csc_form& in_mat) {
    if(this == &in_mat) return *this;
    access::rw(n_rows) = in_mat.n_rows;
    access::rw(n_cols) = in_mat.n_cols;
    access::rw(n_elem) = in_mat.n_elem;
    init(n_elem);
    in_mat.copy_to(row_idx.get(), col_ptr.get(), val_idx.get());
    return *this;
}

template<sp_d data_t, sp_i index_t> csc_form<data_t, index_t>& csc_form<data_t, index_t>::operator=(csc_form&& in_mat) noexcept {
    if(this == &in_mat) return *this;
    access::rw(n_rows) = in_mat.n_rows;
    access::rw(n_cols) = in_mat.n_cols;
    access::rw(n_elem) = in_mat.n_elem;
    row_idx = std::move(in_mat.row_idx);
    col_ptr = std::move(in_mat.col_ptr);
    val_idx = std::move(in_mat.val_idx);
    return *this;
}

template<sp_d data_t, sp_i index_t> void csc_form<data_t, index_t>::print() const {
    suanpan_info("A sparse matrix in triplet form with size of {} by {}, the sparsity of {:.3f}%.\n", n_rows, n_cols, 1E2 - static_cast<double>(n_elem) / static_cast<double>(n_rows) / static_cast<double>(n_cols) * 1E2);
    if(n_elem > index_t(1000)) {
        suanpan_info("More than 1000 elements exist.\n");
        return;
    }

    index_t c_idx = 1;
    for(index_t I = 0; I < n_elem; ++I) {
        if(I >= col_ptr[c_idx]) ++c_idx;
        suanpan_info("({}, {}) ===> {:+.8E}\n", row_idx[I], c_idx - 1, val_idx[I]);
    }
}

template<sp_d data_t, sp_i index_t> template<sp_d in_dt, sp_i in_it> csc_form<data_t, index_t>::csc_form(triplet_form<in_dt, in_it>& in_mat, const SparseBase base, const bool full)
    : n_rows(index_t(in_mat.n_rows))
    , n_cols(index_t(in_mat.n_cols)) {
    if(full) in_mat.full_csc_condense();
    else in_mat.csc_condense();

    init(access::rw(n_elem) = index_t(in_mat.n_elem));

    const sp_i auto shift = index_t(base);

    suanpan::for_each(in_mat.n_elem, [&](const in_it I) {
        row_idx[I] = index_t(in_mat.row_idx[I]) + shift;
        val_idx[I] = data_t(in_mat.val_idx[I]);
    });

    in_it current_pos = 0, current_col = 0;

    while(current_pos < in_mat.n_elem)
        if(in_mat.col_idx[current_pos] < current_col) ++current_pos;
        else col_ptr[current_col++] = index_t(current_pos) + shift;

    col_ptr[0] = shift;
    col_ptr[n_cols] = n_elem + shift;
}

template<sp_d data_t, sp_i index_t> template<sp_d in_dt, sp_i in_it> csc_form<data_t, index_t>& csc_form<data_t, index_t>::operator=(triplet_form<in_dt, in_it>& in_mat) {
    in_mat.csc_condense();

    access::rw(n_rows) = index_t(in_mat.n_rows);
    access::rw(n_cols) = index_t(in_mat.n_cols);

    init(access::rw(n_elem) = index_t(in_mat.n_elem));

    suanpan::for_each(in_mat.n_elem, [&](const in_it I) {
        row_idx[I] = index_t(in_mat.row_idx[I]);
        val_idx[I] = data_t(in_mat.val_idx[I]);
    });

    in_it current_pos = 0, current_col = 0;

    while(current_pos < in_mat.n_elem)
        if(in_mat.col_idx[current_pos] < current_col) ++current_pos;
        else col_ptr[current_col++] = index_t(current_pos);

    col_ptr[0] = index_t(0);
    col_ptr[n_cols] = n_elem;

    return *this;
}

#endif
