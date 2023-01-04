/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#ifdef SUANPAN_MT
#include <tbb/concurrent_set.h>
#include <tbb/concurrent_map.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>

namespace suanpan {
    template<typename T> using vector = tbb::concurrent_vector<T>;
    template<typename T> using set = tbb::concurrent_set<T>;
    template<typename T> using unordered_set = tbb::concurrent_unordered_set<T, std::hash<T>>;
    template<typename T, typename D> using map = tbb::concurrent_map<T, D, std::hash<T>>;
    template<typename T, typename D> using unordered_map = tbb::concurrent_unordered_map<T, D, std::hash<T>>;

    template<typename T> using graph = vector<set<T>>;
}
#else
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>

namespace suanpan {
	template<typename T> using vector = std::vector<T>;
	template<typename T> using set = std::set<T>;
	template<typename T> using unordered_set = std::unordered_set<T>;
	template<typename T, typename D> using map = std::map<T, D>;
	template<typename T, typename D> using unordered_map = std::unordered_map<T, D>;

	template<typename T> using graph = vector<set<T>>;
}
#endif

template<sp_i T> uvec to_uvec(const suanpan::set<T>& in) {
    uvec out(in.size(), fill::none);
    auto I = 0llu;
    for(const auto J : in) out(I++) = static_cast<uword>(J);
    return out;
}

template<sp_i T> uvec to_uvec(const suanpan::unordered_set<T>& in) {
    uvec out(in.size(), fill::none);
    auto I = 0llu;
    for(const auto J : in) out(I++) = static_cast<uword>(J);
    return out;
}

#endif

//! @}
