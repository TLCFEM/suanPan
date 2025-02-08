/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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

#include "sync_ostream.h"
#include <mutex>

namespace suanpan {
    extern std::mutex print_mutex;
}

sync_ostream::sync_ostream(std::ostream& in_ostream)
    : print_lock{suanpan::print_mutex}
    , ostream{&in_ostream} {}

sync_ostream::sync_ostream(sync_ostream&& other) noexcept
    : print_lock{std::move(other.print_lock)}
    , ostream{other.ostream} {}

sync_ostream& sync_ostream::operator=(sync_ostream&& other) noexcept {
    if(this == &other) return *this;
    print_lock = std::move(other.print_lock);
    ostream = other.ostream;
    return *this;
}
