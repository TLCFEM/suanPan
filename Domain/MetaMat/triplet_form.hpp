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

#ifndef TRIPLET_FORM
#define TRIPLET_FORM

#include "comparator.hpp"
#include "sparse_form.hpp"

template<typename data_t, typename index_t> class csc_form;
template<typename data_t, typename index_t> class csr_form;

template<typename data_t, typename index_t> class triplet_form final : public sparse_form<data_t, index_t, triplet_form<data_t, index_t>> {
	using sparse_form<data_t, index_t, triplet_form<data_t, index_t>>::bin;

	void condense() const;

	void copy_memory(index_t, const index_t*, const index_t*, const data_t*) override;
public:
	using sparse_form<data_t, index_t, triplet_form<data_t, index_t>>::n_rows;
	using sparse_form<data_t, index_t, triplet_form<data_t, index_t>>::n_cols;
	using sparse_form<data_t, index_t, triplet_form<data_t, index_t>>::n_elem;
	using sparse_form<data_t, index_t, triplet_form<data_t, index_t>>::c_size;

	index_t* row_idx = nullptr; // index storage
	index_t* col_idx = nullptr; // index storage
	data_t* val_idx = nullptr;  // value storage

	bool csc_sorted = false;
	bool csr_sorted = false;

	triplet_form() = default;
	triplet_form(index_t, index_t, index_t = 0, bool = false, bool = false);

	~triplet_form() override;

	triplet_form(const triplet_form&);                // copy ctor
	triplet_form(triplet_form&&) noexcept;            // move ctor
	triplet_form& operator=(const triplet_form&);     // copy assignment
	triplet_form& operator=(triplet_form&&) noexcept; // move assignment

	template<typename in_dt, typename in_it> triplet_form(const triplet_form<in_dt, in_it>&, int);

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
	void spy() override;

	void csr_condense() const;
	void csc_condense() const;
	[[nodiscard]] bool csr_sort() const;
	[[nodiscard]] bool csc_sort() const;

	[[nodiscard]] index_t row(index_t) const;
	[[nodiscard]] index_t col(index_t) const;
	[[nodiscard]] data_t val(index_t) const;

	const data_t& operator()(index_t, index_t) const;
	data_t& at(index_t, index_t);

	template<typename T2> triplet_form<data_t, index_t> operator*(T2);
	template<typename T2> triplet_form<data_t, index_t> operator/(T2);
	template<typename T2> triplet_form<data_t, index_t>& operator*=(T2);
	template<typename T2> triplet_form<data_t, index_t>& operator/=(T2);

	triplet_form<data_t, index_t> operator+(const triplet_form<data_t, index_t>&);
	triplet_form<data_t, index_t> operator-(const triplet_form<data_t, index_t>&);
	triplet_form<data_t, index_t>& operator+=(const triplet_form<data_t, index_t>&);
	triplet_form<data_t, index_t>& operator-=(const triplet_form<data_t, index_t>&);

	explicit triplet_form(const csc_form<data_t, index_t>&);
	triplet_form& operator=(const csc_form<data_t, index_t>&);
	triplet_form<data_t, index_t> operator+(const csc_form<data_t, index_t>&);
	triplet_form<data_t, index_t> operator-(const csc_form<data_t, index_t>&);
	triplet_form<data_t, index_t>& operator+=(const csc_form<data_t, index_t>&);
	triplet_form<data_t, index_t>& operator-=(const csc_form<data_t, index_t>&);

	explicit triplet_form(const csr_form<data_t, index_t>&);
	triplet_form& operator=(const csr_form<data_t, index_t>&);
	triplet_form<data_t, index_t> operator+(const csr_form<data_t, index_t>&);
	triplet_form<data_t, index_t> operator-(const csr_form<data_t, index_t>&);
	triplet_form<data_t, index_t>& operator+=(const csr_form<data_t, index_t>&);
	triplet_form<data_t, index_t>& operator-=(const csr_form<data_t, index_t>&);

	Mat<data_t> operator*(const Col<data_t>&) override;
	Mat<data_t> operator*(const Mat<data_t>&) override;
};

template<typename data_t, typename index_t> void triplet_form<data_t, index_t>::condense() const {
	if(c_size < 2) return;

	auto last_row = row_idx[0], last_col = col_idx[0];

	index_t current_pos = 0;
	auto last_sum = 0.;

	index_t max_row = 0, max_col = 0;

	for(index_t I = 0; I < c_size; ++I) {
		if(last_row != row_idx[I] || last_col != col_idx[I]) {
			if(last_sum != 0.) {
				row_idx[current_pos] = last_row;
				col_idx[current_pos] = last_col;
				val_idx[current_pos] = last_sum;
				if(last_row > max_row) max_row = last_row;
				if(last_col > max_col) max_col = last_col;
				++current_pos;
				last_sum = 0.;
			}
			last_row = row_idx[I];
			last_col = col_idx[I];
		}
		last_sum += val_idx[I];
	}

	if(last_sum != 0.) {
		row_idx[current_pos] = last_row;
		col_idx[current_pos] = last_col;
		val_idx[current_pos] = last_sum;
		if(last_row > max_row) max_row = last_row;
		if(last_col > max_col) max_col = last_col;
		++current_pos;
	}

	access::rw(n_rows) = max_row + 1;
	access::rw(n_cols) = max_col + 1;
	access::rw(c_size) = current_pos;
}

template<typename data_t, typename index_t> void triplet_form<data_t, index_t>::copy_memory(const index_t size, const index_t* const in_row_idx, const index_t* const in_col_idx, const data_t* const in_val_idx) {
	if(size > n_elem) resize(size);

	std::copy(in_row_idx, in_row_idx + size, this->row_idx);
	std::copy(in_col_idx, in_col_idx + size, this->col_idx);
	std::copy(in_val_idx, in_val_idx + size, this->val_idx);

	access::rw(c_size) = size;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t>::triplet_form(const index_t in_rows, const index_t in_cols, const index_t in_elem, const bool in_csc_sort, const bool in_csr_sort)
	: sparse_form<data_t, index_t, triplet_form<data_t, index_t>>(in_rows, in_cols, in_elem) {
	csc_sorted = in_csc_sort;
	csr_sorted = in_csr_sort;
	triplet_form<data_t, index_t>::init();
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t>::~triplet_form() { triplet_form<data_t, index_t>::reset(); }

template<typename data_t, typename index_t> triplet_form<data_t, index_t>::triplet_form(const triplet_form& in_mat)
	: sparse_form<data_t, index_t, triplet_form<data_t, index_t>>(in_mat.n_rows, in_mat.n_cols, in_mat.n_elem) {
	csc_sorted = in_mat.csc_sorted;
	csr_sorted = in_mat.csr_sorted;
	triplet_form<data_t, index_t>::init();
	triplet_form<data_t, index_t>::copy_memory(in_mat.c_size, in_mat.row_idx, in_mat.col_idx, in_mat.val_idx);
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t>::triplet_form(triplet_form&& in_mat) noexcept {
	triplet_form<data_t, index_t>::reset();
	access::rw(n_rows) = in_mat.n_rows;
	access::rw(n_cols) = in_mat.n_cols;
	access::rw(n_elem) = in_mat.n_elem;
	access::rw(c_size) = in_mat.c_size;
	csc_sorted = in_mat.csc_sorted;
	csr_sorted = in_mat.csr_sorted;
	row_idx = in_mat.row_idx;
	col_idx = in_mat.col_idx;
	val_idx = in_mat.val_idx;
	access::rw(in_mat.n_rows) = access::rw(in_mat.n_cols) = access::rw(in_mat.n_elem) = access::rw(in_mat.c_size) = 0;
	in_mat.row_idx = in_mat.col_idx = nullptr;
	in_mat.val_idx = nullptr;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator=(const triplet_form& in_mat) {
	if(&in_mat != this) {
		init(in_mat.n_rows, in_mat.n_cols, in_mat.n_elem);
		copy_memory(in_mat.c_size, in_mat.row_idx, in_mat.col_idx, in_mat.val_idx);
	}

	return *this;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator=(triplet_form&& in_mat) noexcept {
	triplet_form<data_t, index_t>::reset();
	access::rw(n_rows) = in_mat.n_rows;
	access::rw(n_cols) = in_mat.n_cols;
	access::rw(n_elem) = in_mat.n_elem;
	access::rw(c_size) = in_mat.c_size;
	csc_sorted = in_mat.csc_sorted;
	csr_sorted = in_mat.csr_sorted;
	row_idx = in_mat.row_idx;
	col_idx = in_mat.col_idx;
	val_idx = in_mat.val_idx;
	access::rw(in_mat.n_rows) = access::rw(in_mat.n_cols) = access::rw(in_mat.n_elem) = access::rw(in_mat.c_size) = 0;
	in_mat.row_idx = in_mat.col_idx = nullptr;
	in_mat.val_idx = nullptr;
	return *this;
}

template<typename data_t, typename index_t> template<typename in_dt, typename in_it> triplet_form<data_t, index_t>::triplet_form(const triplet_form<in_dt, in_it>& in_mat, const int base)
	: sparse_form<data_t, index_t, triplet_form<data_t, index_t>>(in_mat.n_rows, in_mat.n_cols, in_mat.n_elem) {
	access::rw(c_size) = in_mat.c_size;
	csc_sorted = in_mat.csc_sorted;
	csr_sorted = in_mat.csr_sorted;

	triplet_form<data_t, index_t>::init();

#ifdef SUANPAN_MT
	tbb::parallel_for(static_cast<index_t>(0), in_mat.c_size, [&](const index_t I) {
		row_idx[I] = index_t(in_mat.row_idx[I]) + base;
		col_idx[I] = index_t(in_mat.col_idx[I]) + base;
		val_idx[I] = data_t(in_mat.val_idx[I]);
	});
#else
	for(index_t I = 0; I < in_mat.c_size; ++I) {
		row_idx[I] = index_t(in_mat.row_idx[I]) + base;
		col_idx[I] = index_t(in_mat.col_idx[I]) + base;
		val_idx[I] = data_t(in_mat.val_idx[I]);
	}
#endif
}

template<typename data_t, typename index_t> const index_t* triplet_form<data_t, index_t>::row_mem() const { return row_idx; }

template<typename data_t, typename index_t> const index_t* triplet_form<data_t, index_t>::col_mem() const { return col_idx; }

template<typename data_t, typename index_t> const data_t* triplet_form<data_t, index_t>::val_mem() const { return val_idx; }

template<typename data_t, typename index_t> void triplet_form<data_t, index_t>::reset() const {
	zeros();
	delete[] row_idx;
	delete[] col_idx;
	delete[] val_idx;
}

template<typename data_t, typename index_t> void triplet_form<data_t, index_t>::zeros() const {
	access::rw(c_size) = 0;
	access::rw(csc_sorted) = false;
	access::rw(csr_sorted) = false;
}

template<typename data_t, typename index_t> data_t triplet_form<data_t, index_t>::max() const { return *std::max_element(val_idx, val_idx + c_size); }

template<typename data_t, typename index_t> bool triplet_form<data_t, index_t>::init() {
	reset();
	if(n_elem == 0) return true;
	row_idx = new(std::nothrow) index_t[n_elem];
	col_idx = new(std::nothrow) index_t[n_elem];
	val_idx = new(std::nothrow) data_t[n_elem];
	if(row_idx == nullptr || col_idx == nullptr || val_idx == nullptr) {
		reset();
		return false;
	}
	return true;
}

template<typename data_t, typename index_t> bool triplet_form<data_t, index_t>::init(const index_t in_elem) {
	if(in_elem <= n_elem) {
		zeros();
		return true;
	}
	access::rw(n_elem) = in_elem;
	return init();
}

template<typename data_t, typename index_t> bool triplet_form<data_t, index_t>::init(const index_t in_row, const index_t in_col, const index_t in_elem) {
	if(n_rows != in_row) access::rw(n_rows) = in_row;
	if(n_cols != in_col) access::rw(n_cols) = in_col;

	return init(in_elem);
}

template<typename data_t, typename index_t> bool triplet_form<data_t, index_t>::resize() {
	const auto copy = *this;

	if(!init(0 == n_elem ? 1 : 2 * n_elem)) return false;

	copy_memory(copy.c_size, copy.row_idx, copy.col_idx, copy.val_idx);

	return true;
}

template<typename data_t, typename index_t> bool triplet_form<data_t, index_t>::resize(const index_t in_elem) {
	const auto copy = *this;

	if(in_elem < c_size || !init(in_elem)) return false;

	copy_memory(copy.c_size, copy.row_idx, copy.col_idx, copy.val_idx);

	return true;
}

template<typename data_t, typename index_t> bool triplet_form<data_t, index_t>::resize(const index_t in_row, const index_t in_col, const index_t in_elem) {
	const auto copy = *this;

	if(in_row < n_rows || in_col < n_cols || in_elem < c_size || !init(in_row, in_col, in_elem)) return false;

	copy_memory(copy.c_size, copy.row_idx, copy.col_idx, copy.val_idx);

	return true;
}

template<typename data_t, typename index_t> void triplet_form<data_t, index_t>::print() const {
	suanpan_info("A sparse matrix in triplet form with size of %u by %u, the density of %.3f.\n", static_cast<unsigned>(n_rows), static_cast<unsigned>(n_cols), static_cast<double>(c_size) / static_cast<double>(n_rows * n_cols) * 100.);
	if(c_size > 1000) {
		suanpan_info("Not going to print all elements as more than 1000 elements exist.\n");
		return;
	}
	for(index_t I = 0; I < c_size; ++I) suanpan_info("(%3u, %3u) ===> %+.4E\n", static_cast<unsigned>(row_idx[I]), static_cast<unsigned>(col_idx[I]), val_idx[I]);
}

template<typename data_t, typename index_t> void triplet_form<data_t, index_t>::spy() {
	if(std::max(n_rows, n_cols) > 100) return;

	csr_condense();

	index_t current_pos = 0;

	for(index_t I = 0; I < n_rows; ++I) {
		for(index_t J = 0; J < n_cols; ++J)
			if(I == row_idx[current_pos] && J == col_idx[current_pos]) {
				suanpan_info(" X");
				++current_pos;
			}
			else suanpan_info(" .");
		suanpan_info("\n");
	}
}

template<typename data_t, typename index_t> void triplet_form<data_t, index_t>::csr_condense() const { if(csr_sort()) condense(); }

template<typename data_t, typename index_t> void triplet_form<data_t, index_t>::csc_condense() const { if(csc_sort()) condense(); }

template<typename data_t, typename index_t> bool triplet_form<data_t, index_t>::csr_sort() const {
	if(c_size < 2) return true;

	const auto new_row_idx = new(std::nothrow) index_t[n_elem];
	const auto new_col_idx = new(std::nothrow) index_t[n_elem];
	const auto new_val_idx = new(std::nothrow) data_t[n_elem];

	if(new_row_idx == nullptr || new_col_idx == nullptr || new_val_idx == nullptr) {
		delete[] new_row_idx;
		delete[] new_col_idx;
		delete[] new_val_idx;
		return false;
	}

	std::vector<index_t> index(c_size);
	for(index_t I = 0; I < static_cast<index_t>(index.size()); ++I) {
		new_row_idx[I] = row_idx[I] * n_cols + col_idx[I];
		index[I] = I;
	}

	suanpan_sort(index.begin(), index.end(), abs_comparator<index_t>(new_row_idx));

#ifdef SUANPAN_MT
	tbb::parallel_for(static_cast<index_t>(0), c_size, [&](const index_t I) {
		new_row_idx[I] = row_idx[index[I]];
		new_col_idx[I] = col_idx[index[I]];
		new_val_idx[I] = val_idx[index[I]];
	});
#else
	for(index_t I = 0; I < c_size; ++I) {
		new_row_idx[I] = row_idx[index[I]];
		new_col_idx[I] = col_idx[index[I]];
		new_val_idx[I] = val_idx[index[I]];
	}
#endif

	reset();

	access::rw(row_idx) = new_row_idx;
	access::rw(col_idx) = new_col_idx;
	access::rw(val_idx) = new_val_idx;
	access::rw(c_size) = static_cast<index_t>(index.size());

	access::rw(csr_sorted) = true;

	return true;
}

template<typename data_t, typename index_t> bool triplet_form<data_t, index_t>::csc_sort() const {
	if(c_size < 2) return true;

	const auto new_row_idx = new(std::nothrow) index_t[n_elem];
	const auto new_col_idx = new(std::nothrow) index_t[n_elem];
	const auto new_val_idx = new(std::nothrow) data_t[n_elem];

	if(new_row_idx == nullptr || new_col_idx == nullptr || new_val_idx == nullptr) {
		delete[] new_row_idx;
		delete[] new_col_idx;
		delete[] new_val_idx;
		return false;
	}

	std::vector<index_t> index(c_size);
	for(index_t I = 0; I < index.size(); ++I) {
		new_row_idx[I] = col_idx[I] * n_rows + row_idx[I];
		index[I] = I;
	}

	suanpan_sort(index.begin(), index.end(), abs_comparator<index_t>(new_row_idx));

#ifdef SUANPAN_MT
	tbb::parallel_for(static_cast<index_t>(0), c_size, [&](const index_t I) {
		new_row_idx[I] = row_idx[index[I]];
		new_col_idx[I] = col_idx[index[I]];
		new_val_idx[I] = val_idx[index[I]];
	});
#else
	for(index_t I = 0; I < c_size; ++I) {
		new_row_idx[I] = row_idx[index[I]];
		new_col_idx[I] = col_idx[index[I]];
		new_val_idx[I] = val_idx[index[I]];
	}
#endif

	reset();

	access::rw(row_idx) = new_row_idx;
	access::rw(col_idx) = new_col_idx;
	access::rw(val_idx) = new_val_idx;
	access::rw(c_size) = static_cast<index_t>(index.size());

	access::rw(csc_sorted) = true;

	return true;
}

template<typename data_t, typename index_t> index_t triplet_form<data_t, index_t>::row(const index_t idx) const {
	if(idx < c_size) return row_idx[idx];
	throw invalid_argument("overflows");
}

template<typename data_t, typename index_t> index_t triplet_form<data_t, index_t>::col(const index_t idx) const {
	if(idx < c_size) return col_idx[idx];
	throw invalid_argument("overflows");
}

template<typename data_t, typename index_t> data_t triplet_form<data_t, index_t>::val(const index_t idx) const {
	if(idx < c_size) return val_idx[idx];
	throw invalid_argument("overflows");
}

template<typename data_t, typename index_t> const data_t& triplet_form<data_t, index_t>::operator()(const index_t row, const index_t col) const {
	for(index_t I = 0; I < c_size; ++I) if(row_idx[I] == row && col_idx[I] == col) return val_idx[I];

	access::rw(bin) = 0.;

	return bin;
}

template<typename data_t, typename index_t> data_t& triplet_form<data_t, index_t>::at(const index_t row, const index_t col) {
	if(row >= n_rows || col >= n_cols) return access::rw(bin);

	if(csr_sorted) csr_sorted = false;
	if(csc_sorted) csc_sorted = false;
	if(c_size == n_elem) resize();
	row_idx[c_size] = row;
	col_idx[c_size] = col;
	return val_idx[access::rw(c_size)++];
}

template<typename data_t, typename index_t> template<typename T2> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::operator*(const T2 scalar) {
	triplet_form<data_t, index_t> copy = *this;

	if(suanpan::approx_equal(1., scalar)) return copy;

#ifdef SUANPAN_MT
	if(suanpan::approx_equal(-1., scalar)) tbb::parallel_for(static_cast<index_t>(0), copy.c_size, [&](const index_t I) { copy.val_idx[I] = -copy.val_idx[I]; });
	else if(suanpan::approx_equal(0., scalar)) tbb::parallel_for(static_cast<index_t>(0), copy.c_size, [&](const index_t I) { copy.val_idx[I] = data_t(0.); });
	else tbb::parallel_for(static_cast<index_t>(0), copy.c_size, [&](const index_t I) { copy.val_idx[I] *= data_t(scalar); });
#else
	if(suanpan::approx_equal(-1., scalar)) for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] = -copy.val_idx[I];
	else if(suanpan::approx_equal(0., scalar)) for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] = data_t(0.);
	else for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] *= data_t(scalar);
#endif

	return copy;
}

template<typename data_t, typename index_t> template<typename T2> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::operator/(const T2 scalar) {
	triplet_form<data_t, index_t> copy = *this;

	if(suanpan::approx_equal(1., scalar)) return copy;

#ifdef SUANPAN_MT
	if(suanpan::approx_equal(-1., scalar)) tbb::parallel_for(static_cast<index_t>(0), copy.c_size, [&](const index_t I) { copy.val_idx[I] = -copy.val_idx[I]; });
	else tbb::parallel_for(static_cast<index_t>(0), copy.c_size, [&](const index_t I) { copy.val_idx[I] /= data_t(scalar); });
#else
	if(suanpan::approx_equal(-1., scalar)) for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] = -copy.val_idx[I];
	else for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] /= data_t(scalar);
#endif

	return copy;
}

template<typename data_t, typename index_t> template<typename T2> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator*=(const T2 scalar) {
	if(suanpan::approx_equal(1., scalar)) return *this;

#ifdef SUANPAN_MT
	if(suanpan::approx_equal(-1., scalar)) tbb::parallel_for(static_cast<index_t>(0), c_size, [&](const index_t I) { val_idx[I] = -val_idx[I]; });
	else if(suanpan::approx_equal(0., scalar)) triplet_form<data_t, index_t>::zeros();
	else tbb::parallel_for(static_cast<index_t>(0), c_size, [&](const index_t I) { val_idx[I] *= data_t(scalar); });
#else
	if(suanpan::approx_equal(-1., scalar)) for(index_t I = 0; I < c_size; ++I) val_idx[I] = -val_idx[I];
	else if(suanpan::approx_equal(0., scalar)) triplet_form<data_t, index_t>::zeros();
	else for(index_t I = 0; I < c_size; ++I) val_idx[I] *= data_t(scalar);
#endif

	return *this;
}

template<typename data_t, typename index_t> template<typename T2> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator/=(const T2 scalar) {
	if(suanpan::approx_equal(1., scalar)) return *this;

#ifdef SUANPAN_MT
	if(suanpan::approx_equal(-1., scalar)) tbb::parallel_for(static_cast<index_t>(0), c_size, [&](const index_t I) { val_idx[I] = -val_idx[I]; });
	else tbb::parallel_for(static_cast<index_t>(0), c_size, [&](const index_t I) { val_idx[I] /= data_t(scalar); });
#else
	if(suanpan::approx_equal(-1., scalar)) for(auto I = 0; I < c_size; ++I) val_idx[I] = -val_idx[I];
	else for(auto I = 0; I < c_size; ++I) val_idx[I] /= data_t(scalar);
#endif

	return *this;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::operator+(const triplet_form<data_t, index_t>& in_mat) {
	triplet_form<data_t, index_t> copy = *this;
	return copy += in_mat;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::operator-(const triplet_form<data_t, index_t>& in_mat) {
	triplet_form<data_t, index_t> copy = *this;
	return copy -= in_mat;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator+=(const triplet_form<data_t, index_t>& in_mat) {
	if(const auto new_size = c_size + in_mat.c_size; n_elem < new_size) resize(new_size);

	for(index_t I = 0; I < in_mat.c_size; ++I) at(in_mat.row(I), in_mat.col(I)) = in_mat.val(I);

	if(in_mat.n_rows > n_rows) access::rw(n_rows) = in_mat.n_rows;
	if(in_mat.n_cols > n_cols) access::rw(n_cols) = in_mat.n_cols;

	return *this;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator-=(const triplet_form<data_t, index_t>& in_mat) {
	if(const auto new_size = c_size + in_mat.c_size; n_elem < new_size) resize(new_size);

	for(index_t I = 0; I < in_mat.c_size; ++I) at(in_mat.row(I), in_mat.col(I)) = -in_mat.val(I);

	if(in_mat.n_rows > n_rows) access::rw(n_rows) = in_mat.n_rows;
	if(in_mat.n_cols > n_cols) access::rw(n_cols) = in_mat.n_cols;

	return *this;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t>::triplet_form(const csc_form<data_t, index_t>& in_mat) { *this = in_mat; }

template<typename data_t, typename index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator=(const csc_form<data_t, index_t>& in_mat) {
	init(in_mat.n_rows, in_mat.n_cols, in_mat.c_size);

	access::rw(c_size) = in_mat.c_size;

	index_t c_idx = 1;
	for(index_t I = 0; I < in_mat.c_size; ++I) {
		if(I >= in_mat.col_ptr[c_idx]) ++c_idx;
		row_idx[I] = in_mat.row_idx[I];
		col_idx[I] = c_idx - 1;
		val_idx[I] = in_mat.val_idx[I];
	}

	return *this;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::operator+(const csc_form<data_t, index_t>& in_mat) {
	triplet_form<data_t, index_t> copy(in_mat);
	return copy += *this;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::operator-(const csc_form<data_t, index_t>& in_mat) {
	triplet_form<data_t, index_t> copy(in_mat);
	return copy -= *this;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator+=(const csc_form<data_t, index_t>& in_mat) { return *this += triplet_form<data_t, index_t>(in_mat); }

template<typename data_t, typename index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator-=(const csc_form<data_t, index_t>& in_mat) { return *this -= triplet_form<data_t, index_t>(in_mat); }

template<typename data_t, typename index_t> triplet_form<data_t, index_t>::triplet_form(const csr_form<data_t, index_t>& in_mat) { *this = in_mat; }

template<typename data_t, typename index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator=(const csr_form<data_t, index_t>& in_mat) {
	init(in_mat.n_rows, in_mat.n_cols, in_mat.c_size);

	access::rw(c_size) = in_mat.c_size;

	index_t c_idx = 1;
	for(index_t I = 0; I < in_mat.c_size; ++I) {
		if(I >= in_mat.row_ptr[c_idx]) ++c_idx;
		row_idx[I] = c_idx - 1;
		col_idx[I] = in_mat.col_idx[I];
		val_idx[I] = in_mat.val_idx[I];
	}

	return *this;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::operator+(const csr_form<data_t, index_t>& in_mat) {
	triplet_form<data_t, index_t> copy(in_mat);
	return copy += *this;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t> triplet_form<data_t, index_t>::operator-(const csr_form<data_t, index_t>& in_mat) {
	triplet_form<data_t, index_t> copy(in_mat);
	return copy -= *this;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator+=(const csr_form<data_t, index_t>& in_mat) { return *this += triplet_form<data_t, index_t>(in_mat); }

template<typename data_t, typename index_t> triplet_form<data_t, index_t>& triplet_form<data_t, index_t>::operator-=(const csr_form<data_t, index_t>& in_mat) { return *this -= triplet_form<data_t, index_t>(in_mat); }

template<typename data_t, typename index_t> Mat<data_t> triplet_form<data_t, index_t>::operator*(const Col<data_t>& in_mat) {
	Mat<data_t> out_mat(in_mat.n_rows, in_mat.n_cols, fill::zeros);

	for(index_t I = 0; I < c_size; ++I) out_mat(row_idx[I]) += val_idx[I] * in_mat(col_idx[I]);

	return out_mat;
}

template<typename data_t, typename index_t> Mat<data_t> triplet_form<data_t, index_t>::operator*(const Mat<data_t>& in_mat) {
	Mat<data_t> out_mat(in_mat.n_rows, in_mat.n_cols, fill::zeros);

	for(index_t I = 0; I < c_size; ++I) out_mat.row(row_idx[I]) += val_idx[I] * in_mat.row(col_idx[I]);

	return out_mat;
}

template<typename data_t, typename index_t> triplet_form<data_t, index_t> operator+(const triplet_form<data_t, index_t>& mat_a, const triplet_form<data_t, index_t>& mat_b) {
	auto out = mat_a;
	out += mat_b;
	return out;
}

#endif
