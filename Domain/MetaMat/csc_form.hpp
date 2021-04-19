/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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

#include "sparse_form.hpp"

template<typename data_t, typename index_t> class triplet_form;
template<typename data_t, typename index_t> class csr_form;

template<typename data_t, typename index_t> class csc_form final : public sparse_form<data_t, index_t, csc_form<data_t, index_t>> {
	using sparse_form<data_t, index_t, csc_form<data_t, index_t>>::bin;

	void copy_memory(index_t, const index_t*, const index_t*, const data_t*) override;
public:
	using sparse_form<data_t, index_t, csc_form<data_t, index_t>>::n_rows;
	using sparse_form<data_t, index_t, csc_form<data_t, index_t>>::n_cols;
	using sparse_form<data_t, index_t, csc_form<data_t, index_t>>::n_elem;
	using sparse_form<data_t, index_t, csc_form<data_t, index_t>>::n_alloc;

	index_t* row_idx = nullptr; // index storage
	index_t* col_ptr = nullptr; // index storage
	data_t* val_idx = nullptr;  // value storage

	csc_form() = default;

	~csc_form() override;

	csc_form(const csc_form&);                // copy ctor
	csc_form(csc_form&&) noexcept;            // move ctor
	csc_form& operator=(const csc_form&);     // copy assignment
	csc_form& operator=(csc_form&&) noexcept; // move assignment

	[[nodiscard]] const index_t* row_mem() const override;
	[[nodiscard]] const index_t* col_mem() const override;
	[[nodiscard]] const data_t* val_mem() const override;

	void reset() const override;
	void zeros() const override;

	[[nodiscard]] data_t max() const override;

	bool init() override;
	bool init(index_t) override;
	bool init(index_t, index_t, index_t) override;
	bool resize() override;
	bool resize(index_t) override;
	bool resize(index_t, index_t, index_t) override;

	void print() const override;

	template<typename T2> csc_form<data_t, index_t> operator*(T2);
	template<typename T2> csc_form<data_t, index_t> operator/(T2);
	template<typename T2> csc_form<data_t, index_t>& operator*=(T2);
	template<typename T2> csc_form<data_t, index_t>& operator/=(T2);

	template<typename in_dt, typename in_it> explicit csc_form(triplet_form<in_dt, in_it>&);
	template<typename in_dt, typename in_it> explicit csc_form(triplet_form<in_dt, in_it>&&);
	template<typename in_dt, typename in_it> explicit csc_form(triplet_form<in_dt, in_it>&, int);
	template<typename in_dt, typename in_it> csc_form& operator=(const triplet_form<in_dt, in_it>&);

	explicit csc_form(const csr_form<data_t, index_t>&);
	csc_form& operator=(const csr_form<data_t, index_t>&);

	const data_t& operator()(index_t, index_t) const;

	Mat<data_t> operator*(const Col<data_t>&) override;
	Mat<data_t> operator*(const Mat<data_t>&) override;
};

template<typename data_t, typename index_t> void csc_form<data_t, index_t>::copy_memory(const index_t in_size, const index_t* const in_row_idx, const index_t* const in_col_ptr, const data_t* const in_val_idx) {
	if(in_size > n_alloc) resize(in_size);

	std::copy(in_row_idx, in_row_idx + in_size, this->row_idx);
	std::copy(in_col_ptr, in_col_ptr + n_cols + 1, this->col_ptr);
	std::copy(in_val_idx, in_val_idx + in_size, this->val_idx);

	access::rw(n_elem) = in_size;
}

template<typename data_t, typename index_t> csc_form<data_t, index_t>::~csc_form() { csc_form<data_t, index_t>::reset(); }

template<typename data_t, typename index_t> csc_form<data_t, index_t>::csc_form(const csc_form& in_mat)
	: sparse_form<data_t, index_t, csc_form<data_t, index_t>>(in_mat.n_rows, in_mat.n_cols, in_mat.n_alloc) {
	csc_form<data_t, index_t>::init();
	csc_form<data_t, index_t>::copy_memory(in_mat.n_elem, in_mat.row_idx, in_mat.col_ptr, in_mat.val_idx);
}

template<typename data_t, typename index_t> csc_form<data_t, index_t>::csc_form(csc_form&& in_mat) noexcept {
	csc_form<data_t, index_t>::reset();
	access::rw(n_rows) = in_mat.n_rows;
	access::rw(n_cols) = in_mat.n_cols;
	access::rw(n_elem) = in_mat.n_elem;
	access::rw(n_alloc) = in_mat.n_alloc;
	row_idx = in_mat.row_idx;
	col_ptr = in_mat.col_ptr;
	val_idx = in_mat.val_idx;
	access::rw(in_mat.n_rows) = access::rw(in_mat.n_cols) = access::rw(in_mat.n_alloc) = access::rw(in_mat.n_elem) = 0;
	in_mat.row_idx = nullptr;
	in_mat.col_ptr = nullptr;
	in_mat.val_idx = nullptr;
}

template<typename data_t, typename index_t> csc_form<data_t, index_t>& csc_form<data_t, index_t>::operator=(const csc_form& in_mat) {
	if(this == &in_mat) return *this;

	init(in_mat.n_rows, in_mat.n_cols, in_mat.n_alloc);
	copy_memory(in_mat.n_elem, in_mat.row_idx, in_mat.col_ptr, in_mat.val_idx);

	return *this;
}

template<typename data_t, typename index_t> csc_form<data_t, index_t>& csc_form<data_t, index_t>::operator=(csc_form&& in_mat) noexcept {
	if(this == &in_mat) return *this;

	reset();
	access::rw(n_rows) = in_mat.n_rows;
	access::rw(n_cols) = in_mat.n_cols;
	access::rw(n_elem) = in_mat.n_elem;
	access::rw(n_alloc) = in_mat.n_alloc;
	row_idx = in_mat.row_idx;
	col_ptr = in_mat.col_ptr;
	val_idx = in_mat.val_idx;
	access::rw(in_mat.n_rows) = access::rw(in_mat.n_cols) = access::rw(in_mat.n_alloc) = access::rw(in_mat.n_elem) = 0;
	in_mat.row_idx = nullptr;
	in_mat.col_ptr = nullptr;
	in_mat.val_idx = nullptr;

	return *this;
}

template<typename data_t, typename index_t> const index_t* csc_form<data_t, index_t>::row_mem() const { return row_idx; }

template<typename data_t, typename index_t> const index_t* csc_form<data_t, index_t>::col_mem() const { return col_ptr; }

template<typename data_t, typename index_t> const data_t* csc_form<data_t, index_t>::val_mem() const { return val_idx; }

template<typename data_t, typename index_t> void csc_form<data_t, index_t>::reset() const {
	zeros();
	delete[] row_idx;
	delete[] col_ptr;
	delete[] val_idx;
}

template<typename data_t, typename index_t> void csc_form<data_t, index_t>::zeros() const { access::rw(n_elem) = 0; }

template<typename data_t, typename index_t> data_t csc_form<data_t, index_t>::max() const { return *std::max_element(val_idx, val_idx + n_elem); }

template<typename data_t, typename index_t> bool csc_form<data_t, index_t>::init() {
	reset();
	row_idx = new(std::nothrow) index_t[n_alloc];
	col_ptr = new(std::nothrow) index_t[n_cols + index_t(1)];
	val_idx = new(std::nothrow) data_t[n_alloc];

	if(row_idx == nullptr || col_ptr == nullptr || val_idx == nullptr) {
		reset();
		return false;
	}
	return true;
}

template<typename data_t, typename index_t> bool csc_form<data_t, index_t>::init(const index_t in_elem) {
	if(in_elem <= n_alloc) {
		zeros();
		return true;
	}
	access::rw(n_alloc) = in_elem;
	return init();
}

template<typename data_t, typename index_t> bool csc_form<data_t, index_t>::init(const index_t in_row, const index_t in_col, const index_t in_elem) {
	if(n_rows != in_row) access::rw(n_rows) = in_row;
	if(n_cols != in_col) access::rw(n_cols) = in_col;

	return init(in_elem);
}

template<typename data_t, typename index_t> bool csc_form<data_t, index_t>::resize() {
	const auto copy = *this;

	if(!init(n_alloc == 0 ? 1 : 2 * n_alloc)) return false;

	copy_memory(copy.n_elem, copy.row_idx, copy.col_ptr, copy.val_idx);

	return true;
}

template<typename data_t, typename index_t> bool csc_form<data_t, index_t>::resize(const index_t in_elem) {
	const auto copy = *this;

	if(in_elem <= n_elem || !init(in_elem)) return false;

	copy_memory(copy.n_elem, copy.row_idx, copy.col_ptr, copy.val_idx);

	return true;
}

template<typename data_t, typename index_t> bool csc_form<data_t, index_t>::resize(const index_t in_row, const index_t in_col, const index_t in_elem) {
	const auto copy = *this;

	if(in_row < n_rows || in_col < n_cols || in_elem < n_elem || !init(in_row, in_col, in_elem)) return false;

	copy_memory(copy.n_elem, copy.row_idx, copy.col_ptr, copy.val_idx);

	return true;
}

template<typename data_t, typename index_t> void csc_form<data_t, index_t>::print() const {
	suanpan_info("A sparse matrix in triplet form with size of %u by %u, the sparsity of %.3f.\n", static_cast<unsigned>(n_rows), static_cast<unsigned>(n_cols), 100. - static_cast<double>(n_elem) / static_cast<double>(n_rows) / static_cast<double>(n_cols) * 100.);
	if(n_elem > 1000) {
		suanpan_info("more than 1000 elements exist.\n");
		return;
	}

	index_t c_idx = 1;
	for(index_t I = 0; I < n_elem; ++I) {
		if(I >= col_ptr[c_idx]) ++c_idx;
		suanpan_info("(%3u, %3u) ===> %+.4E\n", static_cast<unsigned>(row_idx[I]), static_cast<unsigned>(c_idx) - 1, val_idx[I]);
	}
}

template<typename data_t, typename index_t> template<typename T2> csc_form<data_t, index_t> csc_form<data_t, index_t>::operator*(const T2 scalar) {
	csc_form<data_t, index_t> copy = *this;

#ifdef SUANPAN_MT
	tbb::parallel_for(static_cast<index_t>(0), copy.n_elem, [&](const index_t I) { copy.val_idx[I] *= data_t(scalar); });
#else
	for(auto I = 0; I < copy.n_elem; ++I) copy.val_idx[I] *= data_t(scalar);
#endif

	return copy;
}

template<typename data_t, typename index_t> template<typename T2> csc_form<data_t, index_t> csc_form<data_t, index_t>::operator/(const T2 scalar) {
	csc_form<data_t, index_t> copy = *this;

#ifdef SUANPAN_MT
	tbb::parallel_for(static_cast<index_t>(0), copy.n_elem, [&](const index_t I) { copy.val_idx[I] /= data_t(scalar); });
#else
	for(auto I = 0; I < copy.n_elem; ++I) copy.val_idx[I] /= data_t(scalar);
#endif

	return copy;
}

template<typename data_t, typename index_t> template<typename T2> csc_form<data_t, index_t>& csc_form<data_t, index_t>::operator*=(const T2 scalar) {
#ifdef SUANPAN_MT
	tbb::parallel_for(static_cast<index_t>(0), n_elem, [&](const index_t I) { val_idx[I] *= data_t(scalar); });
#else
	for(auto I = 0; I < n_elem; ++I) val_idx[I] *= data_t(scalar);
#endif

	return *this;
}

template<typename data_t, typename index_t> template<typename T2> csc_form<data_t, index_t>& csc_form<data_t, index_t>::operator/=(const T2 scalar) {
#ifdef SUANPAN_MT
	tbb::parallel_for(static_cast<index_t>(0), n_elem, [&](const index_t I) { val_idx[I] /= data_t(scalar); });
#else
	for(auto I = 0; I < n_elem; ++I) val_idx[I] /= data_t(scalar);
#endif

	return *this;
}

template<typename data_t, typename index_t> template<typename in_dt, typename in_it> csc_form<data_t, index_t>::csc_form(triplet_form<in_dt, in_it>& old_mat) { *this = old_mat; }

template<typename data_t, typename index_t> template<typename in_dt, typename in_it> csc_form<data_t, index_t>::csc_form(triplet_form<in_dt, in_it>&& old_mat) { *this = std::forward<triplet_form<in_dt, in_it>>(old_mat); }

template<typename data_t, typename index_t> template<typename in_dt, typename in_it> csc_form<data_t, index_t>::csc_form(triplet_form<in_dt, in_it>& in_mat, const int base) {
	in_mat.csc_condense();

	init(index_t(in_mat.n_rows), index_t(in_mat.n_cols), index_t(in_mat.n_elem));

	if(0 == in_mat.n_elem) return;

	access::rw(n_elem) = index_t(in_mat.n_elem);

	std::transform(in_mat.row_idx, in_mat.row_idx + in_mat.n_elem, row_idx, [&](const in_it I) { return index_t(I) + base; });
	std::transform(in_mat.val_idx, in_mat.val_idx + in_mat.n_elem, val_idx, [](const in_dt I) { return data_t(I); });

	in_it current_pos = 0, current_col = 0;
	while(current_pos < in_mat.n_elem)
		if(in_mat.col_idx[current_pos] < current_col) ++current_pos;
		else col_ptr[current_col++] = index_t(current_pos) + base;

	col_ptr[n_cols] = n_elem + base;
}

template<typename data_t, typename index_t> template<typename in_dt, typename in_it> csc_form<data_t, index_t>& csc_form<data_t, index_t>::operator=(const triplet_form<in_dt, in_it>& in_mat) {
	in_mat.csc_condense();

	init(index_t(in_mat.n_rows), index_t(in_mat.n_cols), index_t(in_mat.n_elem));

	if(0 == in_mat.n_elem) return *this;

	access::rw(n_elem) = index_t(in_mat.n_elem);

	std::transform(in_mat.row_idx, in_mat.row_idx + in_mat.n_elem, row_idx, [](const in_it I) { return index_t(I); });
	std::transform(in_mat.val_idx, in_mat.val_idx + in_mat.n_elem, val_idx, [](const in_dt I) { return data_t(I); });

	in_it current_pos = 0, current_col = 0;
	while(current_pos < in_mat.n_elem)
		if(in_mat.col_idx[current_pos] < current_col) ++current_pos;
		else col_ptr[current_col++] = index_t(current_pos);

	col_ptr[n_cols] = n_elem;

	return *this;
}

template<typename data_t, typename index_t> csc_form<data_t, index_t>::csc_form(const csr_form<data_t, index_t>& in_mat) { *this = in_mat; }

template<typename data_t, typename index_t> csc_form<data_t, index_t>& csc_form<data_t, index_t>::operator=(const csr_form<data_t, index_t>& in_mat) {
	csc_form<data_t, index_t>::init(in_mat.n_rows, in_mat.n_cols, in_mat.n_elem);

#ifdef SUANPAN_MT
	tbb::parallel_for(static_cast<index_t>(0), n_cols + 1, [&](const index_t I) { col_ptr[I] = 0; });
#else
	for(index_t I = 0; I <= n_cols; ++I) col_ptr[I] = 0;
#endif

	for(index_t I = 0; I < in_mat.n_elem; ++I) ++col_ptr[in_mat.col_idx[I] + 1];
	for(index_t I = 2; I <= n_cols; ++I) col_ptr[I] += col_ptr[I - 1];

	std::vector<index_t> counter(n_cols, 0);

	index_t r_idx = 1;
	for(index_t I = 0; I < in_mat.n_elem; ++I) {
		if(I >= in_mat.row_ptr[r_idx]) ++r_idx;
		const auto& c_idx = in_mat.col_idx[I];
		const auto c_pos = counter[c_idx]++ + col_ptr[c_idx];
		row_idx[c_pos] = r_idx - 1;
		val_idx[c_pos] = in_mat.val_idx[I];
	}

	access::rw(n_elem) = in_mat.n_elem;

	return *this;
}

template<typename data_t, typename index_t> const data_t& csc_form<data_t, index_t>::operator()(const index_t in_row, const index_t in_col) const {
	if(in_row < n_rows && in_col < n_cols) for(auto I = col_ptr[in_col]; I < col_ptr[in_col + 1]; ++I) if(row_idx[I] == in_row) return val_idx[I];

	access::rw(bin) = 0.;
	return bin;
}

template<typename data_t, typename index_t> Mat<data_t> csc_form<data_t, index_t>::operator*(const Col<data_t>& in_mat) {
	Mat<data_t> out_mat = arma::zeros<Mat<data_t>>(in_mat.n_rows, 1);

	for(index_t I = 0; I < n_cols; ++I) for(auto J = col_ptr[I]; J < col_ptr[I + 1]; ++J) out_mat(row_idx[J]) += val_idx[J] * in_mat(I);

	return out_mat;
}

template<typename data_t, typename index_t> Mat<data_t> csc_form<data_t, index_t>::operator*(const Mat<data_t>& in_mat) {
	Mat<data_t> out_mat = arma::zeros<Mat<data_t>>(in_mat.n_rows, in_mat.n_cols);

#ifdef SUANPAN_MT
	tbb::parallel_for(static_cast<uword>(0), out_mat.n_cols, [&](const uword K) { for(index_t I = 0; I < n_cols; ++I) for(auto J = col_ptr[I]; J < col_ptr[I + 1]; ++J) out_mat(row_idx[J], K) += val_idx[J] * in_mat(I, K); });
#else
	for(uword K = 0; K < out_mat.n_cols; ++K) for(index_t I = 0; I < n_cols; ++I) for(auto J = col_ptr[I]; J < col_ptr[I + 1]; ++J) out_mat(row_idx[J], K) += val_idx[J] * in_mat(I, K);
#endif

	return out_mat;
}

#endif
