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
/**
 * @class FullMatCUDA
 * @brief A FullMatCUDA class that holds matrices.
 *
 * @author tlc
 * @date 17/04/2021
 * @version 0.1.0
 * @file FullMatCUDA.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef FULLMATCUDA_HPP
#define FULLMATCUDA_HPP

#ifdef SUANPAN_CUDA

template<typename T> class FullMatCUDA final : public FullMat<T> {
public:
	using FullMat<T>::n_cols;
	using FullMat<T>::n_rows;
	using FullMat<T>::memory;
	using FullMat<T>::FullMat;

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;

	int solve_trs(Mat<T>&, const Mat<T>&) override;
};

template<typename T> unique_ptr<MetaMat<T>> FullMatCUDA<T>::make_copy() { return make_unique<FullMatCUDA<T>>(*this); }

template<typename T> int FullMatCUDA<T>::solve(Mat<T>& X, const Mat<T>& B) {
	int INFO = 0;

	return INFO;
}

template<typename T> int FullMatCUDA<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	int INFO = 0;

	return INFO;
}

#endif

#endif

//! @}
