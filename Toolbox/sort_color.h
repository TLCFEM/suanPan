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
 * @fn sort_color
 * @brief A four color sorting algorithm.
 * @author tlc
 * @date 22/02/2020
 * @version 0.1.2
 * @file sort_color.h
 * @addtogroup Utility
 * @{
 */

#ifndef COLOR_H
#define COLOR_H

#include <suanPan.h>

using std::vector;
#ifdef SUANPAN_MT
#include <tbb/concurrent_set.h>
using suanpan_set = tbb::concurrent_set<unsigned, std::less<>>;
#else
#include <set>
using suanpan_set = std::set<unsigned, std::less<>>;
#endif
using suanpan_register = vector<suanpan_set>;

auto sort_color_metis(suanpan_register&, int, char);
vector<vector<unsigned>> sort_color_wp(const suanpan_register&);
vector<vector<unsigned>> sort_color_mis(const suanpan_register&);

#endif

//! @}
