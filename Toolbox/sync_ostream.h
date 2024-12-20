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

#ifndef SYNC_OSTREAM_H
#define SYNC_OSTREAM_H

#include <mutex>
#include <ostream>

class sync_ostream {
    std::unique_lock<std::mutex> print_lock;
    std::ostream* ostream;

public:
    explicit sync_ostream(std::ostream&);
    sync_ostream(const sync_ostream&) = delete;
    sync_ostream(sync_ostream&&) noexcept;
    sync_ostream& operator=(const sync_ostream&) = delete;
    sync_ostream& operator=(sync_ostream&&) noexcept;
    ~sync_ostream() = default;

    template<typename T> sync_ostream& operator<<(T&& item) {
        *ostream << std::forward<T>(item);
        return *this;
    }
};

#endif // SYNC_OSTREAM_H
