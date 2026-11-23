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
#include <atomic>
#include <mutex>
#include <ostream>
#include <string>
#include <suanPan.h>

namespace {
    std::atomic<bool> async_enabled{false};

    thread_pool& get_print_pool() {
        // Constructed on first use, destroyed before std::cout at exit.
        static thread_pool pool(1);
        return pool;
    }
} // namespace

void suanpan::print_message(const std::string& message, std::ostream& stream) {
    // Serialise producers so that the queue order matches the call order.
    const std::scoped_lock lock(print_mutex);

    if(async_enabled.load(std::memory_order_relaxed)) {
        // Exactly one worker consumes the FIFO queue, preserving message order.
        // The pool isolates the worker from throwing tasks.
        get_print_pool().push_task([&stream, message] { stream << message; });
    }
    else {
        stream << message;
    }
}

void suanpan::enable_async_print() {
    const std::scoped_lock lock(print_mutex);
    if(async_enabled.exchange(true)) return;
    get_print_pool();
}

void suanpan::disable_async_print() {
    const std::scoped_lock lock(print_mutex);
    if(!async_enabled.exchange(false)) return;
    get_print_pool().wait_for_tasks();
}

void suanpan::flush_print() {
    const std::scoped_lock lock(print_mutex);
    if(async_enabled.load(std::memory_order_relaxed)) get_print_pool().wait_for_tasks();
}

bool suanpan::is_async_print_enabled() { return async_enabled.load(std::memory_order_relaxed); }

suanpan::AsyncPrintGuard::AsyncPrintGuard() { enable_async_print(); }

suanpan::AsyncPrintGuard::~AsyncPrintGuard() { disable_async_print(); }

//! @}
