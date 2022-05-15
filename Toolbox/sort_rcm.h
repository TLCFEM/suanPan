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
 * @fn sort_rcm
 * @brief A renumber function using RCM algorithm.
 *
 * The function takes both mat and sp_mat.
 *
 * Example Usage:
 *
 * ```cpp
 *     sp_mat test_rcm=sprandn(100000, 100000, 0.00005);
 *     auto R = RCM(test_rcm + test_rcm.t());
 * ```
 *
 * R gives the new numbering order of the original symmetric matrix.
 *
 * @author tlc
 * @date 02/08/2017
 * @version 0.1.2
 * @file sort_rcm.h
 * @addtogroup Utility
 * @{
 */

#ifndef RCM_H
#define RCM_H

#include <Domain/MetaMat/triplet_form.hpp>

using std::vector;

uvec sort_rcm(const vector<uvec>&, const uvec&);

template<typename eT> uvec sort_rcm(const SpMat<eT>&);

template<typename eT> uvec sort_rcm(const Mat<eT>&);

template<typename dt, typename it> uvec sort_rcm(const triplet_form<dt, it>&);

#endif

//! @}
