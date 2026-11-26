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

namespace {
    enum class print_mode : unsigned char { sync,
                                            async,
                                            switching };

    // Guards the mode and serialises every stream write, including worker writes.
    std::mutex mode_mutex;

    std::condition_variable mode_cv;

    // Serialises mode changes so only one drain/flip runs at a time.
    std::mutex switch_mutex;

    print_mode mode = print_mode::sync;

    unsigned depth = 0;

    // Restores the mode and wakes blocked producers, even if a stream throws.
    struct restore_mode {
        std::unique_lock<std::mutex>& lock;
        print_mode next;

        ~restore_mode() {
            mode = next;
            lock.unlock();
            mode_cv.notify_all();
        }
    };

    thread_pool& get_print_pool() {
        static thread_pool pool(1);
        return pool;
    }
} // namespace

namespace suanpan {
    void print_message(const std::string& message, std::ostream& stream) {
        std::unique_lock lock(mode_mutex);
        mode_cv.wait(lock, [] { return print_mode::switching != mode; });

        if(print_mode::async == mode) get_print_pool().push_task([&stream, message] {
            const std::scoped_lock write_lock(mode_mutex);
            // A failing stream must not recurse into the logger.
            try {
                stream << message;
            }
            catch(...) {
            }
        });
        else stream << message;
    }

    void enable_async_print() {
        const std::scoped_lock switch_lock(switch_mutex);

        get_print_pool();

        const std::scoped_lock lock(mode_mutex);
        if(0 == depth++) mode = print_mode::async;
    }

    void disable_async_print() {
        const std::scoped_lock switch_lock(switch_mutex);

        std::unique_lock lock(mode_mutex);
        if(0 == depth) return;

        --depth;
        mode = print_mode::switching;
        lock.unlock();

        get_print_pool().wait_for_tasks();

        lock.lock();
        const restore_mode restore{lock, 0 == depth ? print_mode::sync : print_mode::async};
    }

    void flush_print(std::ostream& stream) {
        const std::scoped_lock switch_lock(switch_mutex);

        std::unique_lock lock(mode_mutex);
        if(print_mode::sync == mode) {
            stream.flush();
            return;
        }

        mode = print_mode::switching;
        lock.unlock();

        get_print_pool().wait_for_tasks();

        lock.lock();
        const restore_mode restore{lock, 0 == depth ? print_mode::sync : print_mode::async};
        stream.flush();
    }

    bool is_async_print_enabled() {
        const std::scoped_lock lock(mode_mutex);
        return depth > 0;
    }
} // namespace suanpan

//! @}
