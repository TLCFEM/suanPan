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
 * @fn eig_solver
 * @brief Select the serial ARPACK or distributed PARPACK eigensolver.
 *
 * @author tlc
 * @date 06/11/2026
 * @version 0.1.0
 * @file eig_solver.h
 * @addtogroup Utility
 * @{
 */

#ifndef EIG_SOLVER_H
#define EIG_SOLVER_H

#ifdef SUANPAN_DISTRIBUTED
#include <Toolbox/parpack.h>
#else
#include <Toolbox/arpack.h>
#endif

#endif

//! @}
