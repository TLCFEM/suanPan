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
 * @file async_print.cpp
 * @brief Implementation of asynchronous screen printing.
 *
 * @author tlc
 * @date 23/11/2026
 * @version 0.1.0
 * @addtogroup Utility
 * @{
 */

#include "async_print.h"

#include <Toolbox/thread_pool.hpp>
#include <suanPan.h>

namespace {
    std::atomic async_enabled{false};
    std::mutex print_mutex;

    thread_pool& get_print_pool() {
        static thread_pool pool(1);
        return pool;
    }
} // namespace

namespace suanpan {
    void print_message(const std::string& message, std::ostream& stream) {
        const std::scoped_lock lock(print_mutex);

        if(async_enabled.load(std::memory_order_relaxed)) get_print_pool().push_task([&stream, message] { stream << message; });
        else stream << message;
    }

    void enable_async_print() {
        const std::scoped_lock lock(print_mutex);
        if(async_enabled.exchange(true)) return;
        get_print_pool();
    }

    void disable_async_print() {
        const std::scoped_lock lock(print_mutex);
        if(!async_enabled.exchange(false)) return;
        get_print_pool().wait_for_tasks();
    }

    void flush_print() {
        const std::scoped_lock lock(print_mutex);
        if(async_enabled.load(std::memory_order_relaxed)) get_print_pool().wait_for_tasks();
    }

    bool is_async_print_enabled() { return async_enabled.load(std::memory_order_relaxed); }
} // namespace suanpan

//! @}
