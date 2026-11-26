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
 * @file async_print.h
 * @brief Asynchronous screen printing support.
 *
 * Printing is synchronous by default, which is required in CLI mode. While
 * asynchronous printing is enabled, a single worker thread writes the queued
 * messages. Enabling, disabling and flushing pause producers and drain the
 * queue before the mode changes, so a write can never overlap a mode change.
 *
 * @author tlc
 * @date 23/11/2026
 * @version 0.1.0
 * @addtogroup Utility
 * @{
 */

#ifndef ASYNC_PRINT_H
#define ASYNC_PRINT_H

#include <iosfwd>

namespace suanpan {
    void enable_async_print();
    void disable_async_print();
    void flush_print(std::ostream& stream);

    [[nodiscard]] bool is_async_print_enabled();

    class AsyncPrintGuard final {
    public:
        AsyncPrintGuard() { enable_async_print(); }
        AsyncPrintGuard(const AsyncPrintGuard&) = delete;
        AsyncPrintGuard(AsyncPrintGuard&&) = delete;
        AsyncPrintGuard& operator=(const AsyncPrintGuard&) = delete;
        AsyncPrintGuard& operator=(AsyncPrintGuard&&) = delete;
        ~AsyncPrintGuard() { disable_async_print(); }
    };
} // namespace suanpan

#endif

//! @}
