/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
 * @fn parpack
 * @brief Provide a wrapper for PARPACK.
 *
 * @author tlc
 * @date 06/11/2026
 * @version 0.1.0
 * @file parpack.h
 * @addtogroup Utility
 * @{
 */

#ifndef PARPACK_H
#define PARPACK_H

#include <Domain/MetaMat/MetaMat.hpp>
#include <memory>

int eig_solve(vec&, mat&, const std::shared_ptr<MetaMat<double>>&, const std::shared_ptr<MetaMat<double>>&, unsigned, const char* = "SM");

int eig_solve(cx_vec&, cx_mat&, const std::shared_ptr<MetaMat<double>>&, const std::shared_ptr<MetaMat<double>>&, unsigned, const char* = "LM");

int eig_solve(cx_vec&, const std::shared_ptr<MetaMat<double>>&, unsigned);

#endif

//! @}
