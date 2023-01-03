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

#include "print.h"
#include <suanPan.h>

namespace suanpan {
    std::mutex print_mutex;
}

void suanpan_info(const char* M, ...) {
    if(!SUANPAN_PRINT) return;

    va_list list;
    va_start(list, M);
    const size_t len = std::vsnprintf(nullptr, 0, M, list);
    va_end(list);
    std::vector<char> vec(len + 1);
    va_start(list, M);
    std::vsnprintf(vec.data(), len + 1, M, list);
    va_end(list);

    const std::scoped_lock lock(suanpan::print_mutex);
    SUANPAN_SYNC_COUT << vec.data();
}
