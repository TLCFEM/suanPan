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
 * Screen output is expensive. The suite of logging helpers declared in
 * suanPan.h formats a message on the calling thread and hands it to
 * suanpan::print_message(). When asynchronous printing is enabled, the actual
 * write is performed by a single dedicated worker thread instead of blocking
 * the caller. A single worker consumes a FIFO queue, thus messages appear in
 * submission order.
 *
 * Asynchronous printing is intended for non-interactive (batch) runs only. In
 * CLI mode the prompt must be written before the next line is read, so printing
 * stays synchronous there.
 *
 * @author tlc
 * @date 23/11/2026
 * @version 0.1.0
 * @addtogroup Utility
 * @{
 */

#ifndef ASYNC_PRINT_H
#define ASYNC_PRINT_H

namespace suanpan {
    /**
     * @brief Enable asynchronous printing.
     *
     * Subsequent messages are queued for the print worker. This function is
     * idempotent and shares suanpan::print_mutex with the logging helpers so
     * that enqueueing and mode switches cannot interleave.
     */
    void enable_async_print();

    /**
     * @brief Disable asynchronous printing.
     *
     * All messages queued so far are drained before this function returns.
     * Any message submitted afterward is written synchronously again.
     */
    void disable_async_print();

    /**
     * @brief Block until all queued messages have been written.
     *
     * Asynchronous printing stays enabled. This is useful before swapping the
     * stream buffer of the output stream.
     */
    void flush_print();

    /**
     * @brief Check whether asynchronous printing is currently active.
     */
    [[nodiscard]] bool is_async_print_enabled();

    /**
     * @class AsyncPrintGuard
     * @brief RAII helper that enables asynchronous printing for a scope.
     *
     * The destructor flushes pending output and restores synchronous printing,
     * including when the scope is left via an exception.
     */
    class AsyncPrintGuard final {
    public:
        AsyncPrintGuard() { enable_async_print(); }
        AsyncPrintGuard(const AsyncPrintGuard&) = delete;
        AsyncPrintGuard& operator=(const AsyncPrintGuard&) = delete;
        ~AsyncPrintGuard() { disable_async_print(); }
    };
} // namespace suanpan

#endif

//! @}
