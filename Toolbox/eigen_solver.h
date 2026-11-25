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
 * @fn eigen_solver
 * @brief Select the serial or distributed eigensolvers.
 *
 * @author tlc
 * @date 06/11/2026
 * @version 0.1.0
 * @file eigen_solver.h
 * @addtogroup Utility
 * @{
 */

#ifndef EIGEN_SOLVER_H
#define EIGEN_SOLVER_H

#ifdef SUANPAN_DISTRIBUTED
#include <Toolbox/parpack.h>
#include <Toolbox/pfeast.h>
#else
#include <Toolbox/arpack.h>
#include <Toolbox/feast.h>
#endif

#endif

//! @}
