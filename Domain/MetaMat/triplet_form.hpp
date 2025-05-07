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

#ifndef TRIPLET_FORM_HPP
#define TRIPLET_FORM_HPP

#include <Toolbox/utility.h>
#include <numeric>

template<sp_d data_t, sp_i index_t> class csc_form;
template<sp_d data_t, sp_i index_t> class csr_form;

enum class SparseBase : short unsigned {
    ZERO,
    ONE
};

template<sp_i index_t> class csr_comparator {
    const index_t* const row_idx;
    const index_t* const col_idx;

public:
    csr_comparator(const index_t* const in_row_idx, const index_t* const in_col_idx)
        : row_idx(in_row_idx)
        , col_idx(in_col_idx) {}

    bool operator()(const index_t idx_a, const index_t idx_b) const {
        if(row_idx[idx_a] == row_idx[idx_b]) return col_idx[idx_a] < col_idx[idx_b];
        return row_idx[idx_a] < row_idx[idx_b];
    }
};

template<sp_i index_t> class csc_comparator {
    const index_t* const row_idx;
    const index_t* const col_idx;

public:
    csc_comparator(const index_t* const in_row_idx, const index_t* const in_col_idx)
        : row_idx(in_row_idx)
        , col_idx(in_col_idx) {}

    bool operator()(const index_t idx_a, const index_t idx_b) const {
        if(col_idx[idx_a] == col_idx[idx_b]) return row_idx[idx_a] < row_idx[idx_b];
        return col_idx[idx_a] < col_idx[idx_b];
    }
};

template<sp_d data_t, sp_i index_t> class triplet_form final {
    const data_t bin = data_t(0);

    using index_ptr = std::unique_ptr<index_t[]>;
    using data_ptr = std::unique_ptr<data_t[]>;

    index_ptr row_idx = nullptr; // index storage
    index_ptr col_idx = nullptr; // index storage
    data_ptr val_idx = nullptr;  // value storage

    bool csc_sorted = false;
    bool csr_sorted = false;

    template<sp_d in_dt, sp_i in_it> void copy_to(const std::unique_ptr<in_it[]>& new_row_idx, const std::unique_ptr<in_it[]>& new_col_idx, const std::unique_ptr<in_dt[]>& new_val_idx, const index_t begin, const index_t row_offset, const index_t col_offset, const data_t scalar) const { copy_to(new_row_idx.get(), new_col_idx.get(), new_val_idx.get(), begin, row_offset, col_offset, scalar); }

    template<sp_d in_dt, sp_i in_it> void copy_to(in_it* const new_row_idx, in_it* const new_col_idx, in_dt* const new_val_idx, const index_t begin, const index_t row_offset, const index_t col_offset, const data_t scalar) const {
        suanpan::for_each(n_elem, [&](const index_t I) {
            new_row_idx[I + begin] = in_it(row_idx[I] + row_offset);
            new_col_idx[I + begin] = in_it(col_idx[I] + col_offset);
            new_val_idx[I + begin] = in_dt(scalar * val_idx[I]);
        });
    }

    /**
     * \brief Ensure the container have the capacity for at least in_elem.
     */
    void reserve(const index_t in_elem) {
        if(in_elem <= n_alloc) return;

        access::rw(n_alloc) = index_t(std::pow(2., std::ceil(std::log2(in_elem)) + 1));

        index_ptr new_row_idx(new index_t[n_alloc]);
        index_ptr new_col_idx(new index_t[n_alloc]);
        data_ptr new_val_idx(new data_t[n_alloc]);

        copy_to(new_row_idx, new_col_idx, new_val_idx, 0, 0, 0, 1);

        row_idx = std::move(new_row_idx);
        col_idx = std::move(new_col_idx);
        val_idx = std::move(new_val_idx);
    }

    void invalidate_sorting_flag() { csc_sorted = csr_sorted = false; }

    void condense(bool = false);

    void populate_diagonal() {
        const auto t_elem = std::min(n_rows, n_cols);
        reserve(n_elem + t_elem);
        suanpan::for_each(t_elem, [&](const index_t I) {
            row_idx[n_elem + I] = I;
            col_idx[n_elem + I] = I;
            val_idx[n_elem + I] = data_t(0);
        });
        access::rw(n_elem) += t_elem;
        invalidate_sorting_flag();
    }

public:
    using data_type = data_t;
    using index_type = index_t;

    template<sp_d in_dt, sp_i in_it> friend class csc_form;
    template<sp_d in_dt, sp_i in_it> friend class csr_form;
    template<sp_d in_dt, sp_i in_it> friend class triplet_form;

    const index_t n_rows = 0;
    const index_t n_cols = 0;
    const index_t n_elem = 0;
    const index_t n_alloc = 0;

    triplet_form() = default;
    triplet_form(const triplet_form&);
    triplet_form(triplet_form&&) noexcept;
    triplet_form& operator=(const triplet_form&);
    triplet_form& operator=(triplet_form&&) noexcept;
    ~triplet_form() = default;

    triplet_form(const index_t in_rows, const index_t in_cols, const index_t in_elem = index_t(0))
        : n_rows(in_rows)
        , n_cols(in_cols) { init(in_elem); }

    template<sp_d in_dt> explicit triplet_form(const SpMat<in_dt>&);
    template<sp_d in_dt, sp_i in_it> explicit triplet_form(triplet_form<in_dt, in_it>&, SparseBase = SparseBase::ZERO, bool = false);

    template<sp_d in_dt, sp_i in_it> explicit triplet_form(triplet_form<in_dt, in_it>&& in_mat, const SparseBase in_base = SparseBase::ZERO, const bool in_full = false)
        : triplet_form(in_mat, in_base, in_full) {}

    [[nodiscard]] const index_t* row_mem() const { return row_idx.get(); }

    [[nodiscard]] const index_t* col_mem() const { return col_idx.get(); }

    [[nodiscard]] const data_t* val_mem() const { return val_idx.get(); }

    [[nodiscard]] index_t* row_mem() { return row_idx.get(); }

    [[nodiscard]] index_t* col_mem() { return col_idx.get(); }

    [[nodiscard]] data_t* val_mem() { return val_idx.get(); }

    [[nodiscard]] index_t row(const index_t I) const { return row_idx[I]; }

    [[nodiscard]] index_t col(const index_t I) const { return col_idx[I]; }

    [[nodiscard]] data_t val(const index_t I) const { return val_idx[I]; }

    [[nodiscard]] bool is_csr_sorted() const { return csr_sorted; }

    [[nodiscard]] bool is_csc_sorted() const { return csc_sorted; }

    [[nodiscard]] bool is_empty() const { return 0 == n_elem; }

    [[nodiscard]] data_t max() const {
        if(is_empty()) return data_t(0);
        return *suanpan::max_element(val_idx.get(), val_idx.get() + n_elem);
    }

    void zeros() {
        access::rw(n_elem) = 0;
        invalidate_sorting_flag();
    }

    void init(const index_t in_elem) {
        zeros();
        reserve(in_elem);
    }

    void init(const index_t in_rows, const index_t in_cols, const index_t in_elem) {
        access::rw(n_rows) = in_rows;
        access::rw(n_cols) = in_cols;
        init(in_elem);
    }

    data_t operator()(const index_t row, const index_t col) const {
        for(index_t I = 0; I < n_elem; ++I)
            if(row == row_idx[I] && col == col_idx[I]) return val_idx[I];
        return access::rw(bin) = 0.;
    }

    data_t& at(index_t, index_t);

    auto nullify(const index_t idx) {
        suanpan::for_each(n_elem, [&](const auto I) {
            if(idx == row(I) || idx == col(I)) val_mem()[I] = data_t(0);
        });
    }

    [[nodiscard]] auto extract_col(const index_t idx) {
        csc_condense();

        SpMat<data_t> output(n_rows, 1);
        for(index_t I{0}; I < n_elem; ++I) {
            if(col(I) > idx) break;
            if(col(I) == idx) output.at(row(I), 0) = val(I);
        }
        return output;
    }

    /**
     * @brief Adjusts the size of the container by reserving space and updating the element count.
     *
     * This function reserves memory for the specified number of elements and updates the internal
     * element count to match the provided value.
     *
     * This function is dangerous and should be used with caution.
     * It does not initialize the new elements, and the caller needs to fill them with valid data.
     * Otherwise, the behavior is undefined.
     *
     * @param in_elem The number of elements to reserve and set as the new size.
     */
    void hack_size(const index_t in_elem) {
        reserve(in_elem);
        access::rw(n_elem) = in_elem;
    }

    void print() const;

    void save(std::string_view) const;

    void csr_sort();
    void csc_sort();

    void csr_condense() {
        csr_sort();
        condense(false);
    }

    void csc_condense() {
        csc_sort();
        condense(false);
    }

    void full_csr_condense() {
        populate_diagonal();
        csr_sort();
        condense(true);
    }

    void full_csc_condense() {
        populate_diagonal();
        csc_sort();
        condense(true);
    }

    void assemble(const Mat<data_t>&, const Col<uword>&);
    template<sp_d in_dt, sp_i in_it> void assemble(const triplet_form<in_dt, in_it>&, index_t, index_t, data_t);

    template<sp_d in_dt, sp_i in_it> void assemble(const triplet_form<in_dt, in_it>& in_mat, const std::vector<index_t>& row_shift, const std::vector<index_t>& col_shift, const std::vector<data_t>& scalar) {
        suanpan_assert([&] { if(scalar.size() != row_shift.size() || scalar.size() != col_shift.size()) throw std::invalid_argument("size mismatch detected"); });

        reserve(n_elem + index_t(scalar.size()) * index_t(in_mat.n_elem));

        for(size_t I = 0; I < scalar.size(); ++I) assemble(in_mat, row_shift[I], col_shift[I], scalar[I]);
    }

    Mat<data_t> operator*(const Col<data_t>& in_mat) const {
        Mat<data_t> out_mat(in_mat.n_rows, in_mat.n_cols, fill::zeros);

        for(index_t I = 0; I < n_elem; ++I) out_mat(row_idx[I]) += val_idx[I] * in_mat(col_idx[I]);

        return out_mat;
    }

    Mat<data_t> operator*(const Mat<data_t>& in_mat) const {
        Mat<data_t> out_mat(in_mat.n_rows, in_mat.n_cols, fill::zeros);

        for(index_t I = 0; I < n_elem; ++I) out_mat.row(row_idx[I]) += val_idx[I] * in_mat.row(col_idx[I]);

        return out_mat;
    }

    template<sp_d T2> triplet_form operator*(T2) const;
    template<sp_d T2> triplet_form operator/(T2) const;
    template<sp_d T2> triplet_form& operator*=(T2);
    template<sp_d T2> triplet_form& operator/=(T2);

    triplet_form operator+(const triplet_form& in_mat) const {
        triplet_form copy = *this;
        // ReSharper disable once CppDFAUnusedValue
        return copy += in_mat;
    }

    triplet_form operator-(const triplet_form& in_mat) const {
        triplet_form copy = *this;
        // ReSharper disable once CppDFAUnusedValue
        return copy -= in_mat;
    }

    triplet_form& operator+=(const triplet_form&);
    triplet_form& operator-=(const triplet_form&);

    [[nodiscard]] Col<data_t> diag() const;
    [[nodiscard]] triplet_form diagonal() const;
    [[nodiscard]] triplet_form strictly_upper() const;
    [[nodiscard]] triplet_form strictly_lower() const;
    [[nodiscard]] triplet_form upper() const;
    [[nodiscard]] triplet_form lower() const;
};

template<sp_d data_t, sp_i index_t> void triplet_form<data_t, index_t>::condense(const bool full) {
    if(n_elem < 2) return;

    auto last_row = row_idx[0], last_col = col_idx[0];

    sp_i auto current_pos = index_t(0);
    sp_d auto last_sum = data_t(0);

    auto populate = [&] {
        if(suanpan::approx_equal(last_sum, data_t(0)) && (!full || last_row != last_col)) return;
        row_idx[current_pos] = last_row;
        col_idx[current_pos] = last_col;
        val_idx[current_pos] = last_sum;
        ++current_pos;
        last_sum = data_t(0);
    };

    for(index_t I = 0; I < n_elem; ++I) {
        if(last_row != row_idx[I] || last_col != col_idx[I]) {
            populate();
            last_row = row_idx[I];
            last_col = col_idx[I];
        }
        last_sum += val_idx[I];
    }

    populate();

    access::rw(n_elem) = current_pos;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t>::triplet_form(const triplet_form& in_mat)
    : csc_sorted{in_mat.csc_sorted}
    , csr_sorted{in_mat.csr_sorted}
    , n_rows{in_mat.n_rows}
    , n_cols{in_mat.n_cols} {
    init(in_mat.n_alloc);
    in_mat.copy_to(row_idx, col_idx, val_idx, 0, 0, 0, 1);
    access::rw(n_elem) = in_mat.n_elem;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t>::triplet_form(triplet_form&& in_mat) noexcept
    : row_idx{std::move(in_mat.row_idx)}
    , col_idx{std::move(in_mat.col_idx)}
    , val_idx{std::move(in_mat.val_idx)}
    , csc_sorted{in_mat.csc_sorted}
    , csr_sorted{in_mat.csr_sorted}
    , n_rows{in_mat.n_rows}
    , n_cols{in_mat.n_cols}
    , n_elem{in_mat.n_elem}
    , n_alloc{in_mat.n_alloc} {}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator=(const triplet_form& in_mat) {
    if(this == &in_mat) return *this;
    csc_sorted = in_mat.csc_sorted;
    csr_sorted = in_mat.csr_sorted;
    access::rw(n_rows) = in_mat.n_rows;
    access::rw(n_cols) = in_mat.n_cols;
    init(in_mat.n_alloc);
    in_mat.copy_to(row_idx, col_idx, val_idx, 0, 0, 0, 1);
    access::rw(n_elem) = in_mat.n_elem;
    return *this;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator=(triplet_form&& in_mat) noexcept {
    if(this == &in_mat) return *this;
    csc_sorted = in_mat.csc_sorted;
    csr_sorted = in_mat.csr_sorted;
    access::rw(n_rows) = in_mat.n_rows;
    access::rw(n_cols) = in_mat.n_cols;
    access::rw(n_elem) = in_mat.n_elem;
    access::rw(n_alloc) = in_mat.n_alloc;
    row_idx = std::move(in_mat.row_idx);
    col_idx = std::move(in_mat.col_idx);
    val_idx = std::move(in_mat.val_idx);
    return *this;
}

template<sp_d data_t, sp_i index_t> template<sp_d in_dt> triplet_form<data_t, index_t>::triplet_form(const SpMat<in_dt>& in_mat)
    : n_rows(index_t(in_mat.n_rows))
    , n_cols(index_t(in_mat.n_cols)) {
    init(index_t(in_mat.n_nonzero));

    for(auto I = in_mat.begin(); I != in_mat.end(); ++I) at(I.row(), I.col()) = *I;
}

template<sp_d data_t, sp_i index_t> template<sp_d in_dt, sp_i in_it> triplet_form<data_t, index_t>::triplet_form(triplet_form<in_dt, in_it>& in_mat, const SparseBase base, const bool full)
    : csc_sorted(in_mat.csc_sorted)
    , csr_sorted(in_mat.csr_sorted)
    , n_rows(index_t(in_mat.n_rows))
    , n_cols(index_t(in_mat.n_cols)) {
    if(in_mat.is_empty()) return;

    init(index_t(in_mat.n_alloc));

    if(full) in_mat.full_csc_condense();
    else in_mat.csc_condense();

    const sp_i auto shift = index_t(base);

    in_mat.copy_to(row_idx, col_idx, val_idx, 0, shift, shift, 1);

    access::rw(n_elem) = index_t(in_mat.n_elem);
}

template<sp_d data_t, sp_i index_t> data_t& triplet_form<data_t, index_t>::at(const index_t row, const index_t col) {
    suanpan_assert([&] { if(row >= n_rows || col >= n_cols) throw std::invalid_argument("inconsistent size"); });

    invalidate_sorting_flag();
    reserve(n_elem + 1);
    row_idx[n_elem] = row;
    col_idx[n_elem] = col;
    return val_idx[access::rw(n_elem)++] = data_t(0);
}

template<sp_d data_t, sp_i index_t> void triplet_form<data_t, index_t>::print() const {
    suanpan_info("A sparse matrix in triplet form with size of {} by {}, the density of {:.3f}%.\n", n_rows, n_cols, static_cast<double>(n_elem) / static_cast<double>(n_rows) / static_cast<double>(n_cols) * 1E2);
    if(n_elem > index_t(1000)) {
        suanpan_info("More than 1000 elements exist.\n");
        return;
    }
    for(index_t I = 0; I < n_elem; ++I)
        suanpan_info("({}, {}) ===> {:+.8E}\n", row_idx[I], col_idx[I], val_idx[I]);
}

/**
 * @brief Saves the triplet form matrix to a file.
 *
 * This function writes the matrix data stored in triplet form to a file.
 * The file will contain the number of rows, columns, and non-zero elements
 * on the first line, followed by the row indices, column indices, and values
 * of the non-zero elements on subsequent lines.
 * This follows the Matrix Market (MTX) format.
 *
 * @tparam data_t The data type of the matrix values.
 * @tparam index_t The data type of the row and column indices.
 * @param file_name The name of the file to save the matrix to.
 *
 * @note If the file cannot be opened, the function will return without
 * performing any operation.
 */
template<sp_d data_t, sp_i index_t> void triplet_form<data_t, index_t>::save(const std::string_view file_name) const {
    std::ofstream file(file_name.data());
    if(!file.is_open()) return;
    file << suanpan::format("{} {} {}\n", n_rows, n_cols, n_elem);
    for(index_t I = 0; I < n_elem; ++I) file << suanpan::format("{} {} {}\n", row_idx[I], col_idx[I], val_idx[I]);
    file.close();
}

/**
 * @brief Sorts the COO format into the CSR (Compressed Sparse Row) order.
 *
 * This function ensures that the COO representation of the matrix is sorted
 * by row indices and then by column indices within each row.
 *
 * @tparam data_t The data type of the matrix values.
 * @tparam index_t The data type of the matrix indices.
 *
 * @note After calling this function:
 * - The triplet form will be sorted.
 * - Need to further call `condense()` to remove duplicate entries.
 */
template<sp_d data_t, sp_i index_t> void triplet_form<data_t, index_t>::csr_sort() {
    if(csr_sorted) return;

    std::vector<index_t> index(n_elem);
    std::iota(index.begin(), index.end(), index_t(0));

    suanpan_sort(index.begin(), index.end(), csr_comparator<index_t>(row_idx.get(), col_idx.get()));

    index_ptr new_row_idx(new index_t[n_alloc]);
    index_ptr new_col_idx(new index_t[n_alloc]);
    data_ptr new_val_idx(new data_t[n_alloc]);

    suanpan::for_each(n_elem, [&](const index_t I) {
        new_row_idx[I] = row_idx[index[I]];
        new_col_idx[I] = col_idx[index[I]];
        new_val_idx[I] = val_idx[index[I]];
    });

    row_idx = std::move(new_row_idx);
    col_idx = std::move(new_col_idx);
    val_idx = std::move(new_val_idx);

    csr_sorted = true;
    csc_sorted = false;
}

/**
 * @brief Sorts the COO format into the CSC (Compressed Sparse Column) order.
 *
 * This function ensures that the COO representation of the matrix is sorted
 * by column indices and then by row indices within each column.
 *
 * @tparam data_t The data type of the matrix values.
 * @tparam index_t The data type of the matrix indices.
 *
 * @note After calling this function:
 * - The triplet form will be sorted.
 * - Need to further call `condense()` to remove duplicate entries.
 */
template<sp_d data_t, sp_i index_t> void triplet_form<data_t, index_t>::csc_sort() {
    if(csc_sorted) return;

    std::vector<index_t> index(n_elem);
    std::iota(index.begin(), index.end(), index_t(0));

    suanpan_sort(index.begin(), index.end(), csc_comparator<index_t>(row_idx.get(), col_idx.get()));

    index_ptr new_row_idx(new index_t[n_alloc]);
    index_ptr new_col_idx(new index_t[n_alloc]);
    data_ptr new_val_idx(new data_t[n_alloc]);

    suanpan::for_each(n_elem, [&](const index_t I) {
        new_row_idx[I] = row_idx[index[I]];
        new_col_idx[I] = col_idx[index[I]];
        new_val_idx[I] = val_idx[index[I]];
    });

    row_idx = std::move(new_row_idx);
    col_idx = std::move(new_col_idx);
    val_idx = std::move(new_val_idx);

    csc_sorted = true;
    csr_sorted = false;
}

template<sp_d data_t, sp_i index_t> void triplet_form<data_t, index_t>::assemble(const Mat<data_t>& in_mat, const Col<uword>& in_dof) {
    if(in_mat.empty()) return;

    invalidate_sorting_flag();

    const auto t_elem = n_elem + index_t(in_mat.n_elem);

    reserve(t_elem);

    suanpan::for_each(in_mat.n_elem, [&](const uword I) {
        row_idx[n_elem + I] = index_t(in_dof(I % in_dof.n_elem));
        col_idx[n_elem + I] = index_t(in_dof(I / in_dof.n_elem));
        val_idx[n_elem + I] = in_mat(I);
    });

    access::rw(n_elem) = t_elem;
}

template<sp_d data_t, sp_i index_t> template<sp_d in_dt, sp_i in_it> void triplet_form<data_t, index_t>::assemble(const triplet_form<in_dt, in_it>& in_mat, const index_t row_shift, const index_t col_shift, const data_t scalar) {
    if(in_mat.is_empty()) return;

    invalidate_sorting_flag();

    const auto t_elem = n_elem + in_mat.n_elem;

    reserve(t_elem);

    in_mat.copy_to(row_idx, col_idx, val_idx, n_elem, row_shift, col_shift, scalar);

    access::rw(n_elem) = t_elem;
}

template<sp_d data_t, sp_i index_t> template<sp_d T2> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::operator*(const T2 scalar) const {
    if(suanpan::approx_equal(T2(0), scalar)) return triplet_form(n_rows, n_cols);

    auto copy = *this;

    if(suanpan::approx_equal(T2(1), scalar)) return copy;

    if(suanpan::approx_equal(T2(-1), scalar))
        suanpan_for_each(copy.val_idx.get(), copy.val_idx.get() + copy.n_elem, [](data_t& I) { I = -I; });
    else
        suanpan_for_each(copy.val_idx.get(), copy.val_idx.get() + copy.n_elem, [=](data_t& I) { I *= data_t(scalar); });

    return copy;
}

template<sp_d data_t, sp_i index_t> template<sp_d T2> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::operator/(const T2 scalar) const {
    auto copy = *this;

    if(suanpan::approx_equal(T2(1), scalar)) return copy;

    if(suanpan::approx_equal(T2(-1), scalar))
        suanpan_for_each(copy.val_idx.get(), copy.val_idx.get() + copy.n_elem, [](data_t& I) { I = -I; });
    else
        suanpan_for_each(copy.val_idx.get(), copy.val_idx.get() + copy.n_elem, [=](data_t& I) { I /= data_t(scalar); });

    return copy;
}

template<sp_d data_t, sp_i index_t> template<sp_d T2> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator*=(const T2 scalar) {
    if(suanpan::approx_equal(T2(1), scalar)) return *this;

    if(suanpan::approx_equal(T2(0), scalar)) zeros();
    else if(suanpan::approx_equal(T2(-1), scalar))
        suanpan_for_each(val_idx.get(), val_idx.get() + n_elem, [](data_t& I) { I = -I; });
    else
        suanpan_for_each(val_idx.get(), val_idx.get() + n_elem, [=](data_t& I) { I *= data_t(scalar); });

    return *this;
}

template<sp_d data_t, sp_i index_t> template<sp_d T2> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator/=(const T2 scalar) {
    if(suanpan::approx_equal(T2(1), scalar)) return *this;

    if(suanpan::approx_equal(T2(-1), scalar))
        suanpan_for_each(val_idx.get(), val_idx.get() + n_elem, [](data_t& I) { I = -I; });
    else
        suanpan_for_each(val_idx.get(), val_idx.get() + n_elem, [=](data_t& I) { I /= data_t(scalar); });

    return *this;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator+=(const triplet_form& in_mat) {
    if(in_mat.is_empty()) return *this;

    invalidate_sorting_flag();

    const auto t_elem = n_elem + index_t(in_mat.n_elem);

    reserve(t_elem);

    in_mat.copy_to(row_idx, col_idx, val_idx, n_elem, 0, 0, 1);

    access::rw(n_elem) = t_elem;

    return *this;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator-=(const triplet_form& in_mat) {
    if(in_mat.is_empty()) return *this;

    invalidate_sorting_flag();

    const auto t_elem = n_elem + index_t(in_mat.n_elem);

    reserve(t_elem);

    in_mat.copy_to(row_idx, col_idx, val_idx, n_elem, 0, 0, -1);

    access::rw(n_elem) = t_elem;

    return *this;
}

template<sp_d data_t, sp_i index_t> Col<data_t> triplet_form<data_t, index_t>::diag() const {
    Col<data_t> diag_vec(std::min(n_rows, n_cols), fill::zeros);
    for(index_t I = 0; I < n_elem; ++I)
        if(row(I) == col(I)) diag_vec(row(I)) += val(I);
    return diag_vec;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::diagonal() const {
    auto out_mat = *this;

    suanpan::for_each(out_mat.n_elem, [&](const index_t I) { out_mat.val_idx[I] *= out_mat.row(I) == out_mat.col(I); });

    return out_mat;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::strictly_upper() const {
    auto out_mat = *this;

    suanpan::for_each(out_mat.n_elem, [&](const index_t I) { out_mat.val_idx[I] *= out_mat.row(I) < out_mat.col(I); });

    return out_mat;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::strictly_lower() const {
    auto out_mat = *this;

    suanpan::for_each(out_mat.n_elem, [&](const index_t I) { out_mat.val_idx[I] *= out_mat.col(I) < out_mat.row(I); });

    return out_mat;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::upper() const {
    auto out_mat = *this;

    suanpan::for_each(out_mat.n_elem, [&](const index_t I) { out_mat.val_idx[I] *= out_mat.row(I) <= out_mat.col(I); });

    return out_mat;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::lower() const {
    auto out_mat = *this;

    suanpan::for_each(out_mat.n_elem, [&](const index_t I) { out_mat.val_idx[I] *= out_mat.col(I) <= out_mat.row(I); });

    return out_mat;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> operator+(const triplet_form<data_t, index_t>& mat_a, const triplet_form<data_t, index_t>& mat_b) {
    auto out = mat_a;
    out += mat_b;
    return out;
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> operator+(triplet_form<data_t, index_t>&& mat_a, triplet_form<data_t, index_t>&& mat_b) {
    mat_a += mat_b;
    return std::move(mat_a);
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> operator+(const triplet_form<data_t, index_t>& mat_a, triplet_form<data_t, index_t>&& mat_b) {
    mat_b += mat_a;
    return std::move(mat_b);
}

template<sp_d data_t, sp_i index_t> triplet_form<data_t, index_t> operator+(triplet_form<data_t, index_t>&& mat_a, const triplet_form<data_t, index_t>& mat_b) {
    mat_a += mat_b;
    return std::move(mat_a);
}

#endif
