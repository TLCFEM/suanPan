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
 * @author tlc
 * @date 30/05/2022
 * @version 0.1.0
 * @file container.h
 * @addtogroup Utility
 * @{
 */

#ifndef CONTAINER_H
#define CONTAINER_H

#include <suanPan.h>

using std::vector;

#ifdef SUANPAN_MT
#include <tbb/concurrent_set.h>
#include <tbb/concurrent_unordered_set.h>

namespace suanpan {
    template<typename T> using set = tbb::concurrent_set<T>;
    template<typename T> using unordered_set = tbb::concurrent_unordered_set<T>;
    template<typename T> using graph = vector<set<T>>;
}
#else
#include <set>
#include <unordered_set>

namespace suanpan {
    template<typename T> using set = std::set<T>;
    template<typename T> using unordered_set = std::unordered_set<T>;
    template<typename T> using graph = vector<set<T>>;
}
#endif

#endif

//! @}
