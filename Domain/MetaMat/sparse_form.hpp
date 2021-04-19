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

// ReSharper disable CppUnusedIncludeDirective
#ifndef SPARSE_FORM_HPP
#define SPARSE_FORM_HPP

#include <suanPan.h>
#include <armadillo/armadillo>
#include <vector>
#include <cstring>
#include <Toolbox/utility.h>

template<typename data_t, typename index_t, typename form_t> class sparse_form {
protected:
	const data_t bin = 0.; // bin for out of bound elements

	virtual void copy_memory(index_t, const index_t*, const index_t*, const data_t*) = 0;
public:
	typedef data_t data_type;
	typedef index_t index_type;
	typedef form_t form_type;

	const index_t n_rows = 0;  // number of rows
	const index_t n_cols = 0;  // number of cols
	const index_t n_elem = 0;  // current number of valid elements
	const index_t n_alloc = 0; // maximum number of elements

	sparse_form() = default;
	sparse_form(index_t, index_t, index_t = 0);

	virtual ~sparse_form() = default;

	sparse_form(const sparse_form&) = delete;            // copy ctor
	sparse_form(sparse_form&&) = delete;                 // move ctor
	sparse_form& operator=(const sparse_form&) = delete; // copy assignment
	sparse_form& operator=(sparse_form&&) = delete;      // move assignment

	[[nodiscard]] virtual const index_t* row_mem() const = 0;
	[[nodiscard]] virtual const index_t* col_mem() const = 0;
	[[nodiscard]] virtual const data_t* val_mem() const = 0;

	[[nodiscard]] virtual bool is_empty() const;

	virtual void reset() const = 0;
	virtual void zeros() const = 0;

	[[nodiscard]] virtual data_t max() const = 0;

	virtual bool init() = 0;
	virtual bool init(index_t) = 0;
	virtual bool init(index_t, index_t, index_t) = 0;
	virtual bool resize() = 0;
	virtual bool resize(index_t) = 0;
	virtual bool resize(index_t, index_t, index_t) = 0;

	virtual void print() const;

	virtual Mat<data_t> operator*(const Col<data_t>&) = 0;
	virtual Mat<data_t> operator*(const Mat<data_t>&) = 0;
};

template<typename data_t, typename index_t, typename form_t> sparse_form<data_t, index_t, form_t>::sparse_form(const index_t in_rows, const index_t in_cols, const index_t in_alloc)
	: n_rows(in_rows)
	, n_cols(in_cols)
	, n_alloc(in_alloc) {}

template<typename data_t, typename index_t, typename form_t> bool sparse_form<data_t, index_t, form_t>::is_empty() const { return 0 == n_elem; }

template<typename data_t, typename index_t, typename form_t> void sparse_form<data_t, index_t, form_t>::print() const {}

#endif
