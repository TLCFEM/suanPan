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
 * @file TestAsyncPrint.cpp
 * @brief Tests for the asynchronous screen printer.
 *
 * @author tlc
 * @date 23/11/2026
 * @version 0.1.0
 * @addtogroup Utility
 * @{
 */

#include "CatchHeader.h"

#include <Toolbox/async_print.h>
#include <sstream>
#include <suanPan.h>
#include <thread>

namespace {
    struct counting_buf final : std::stringbuf {
        int sync_count = 0;

        int sync() override {
            ++sync_count;
            return std::stringbuf::sync();
        }
    };

    std::thread::id main_thread_id;
    std::atomic<bool> formatted_on_worker{false};

    struct thread_probe {};
} // namespace

template<> struct fmt::formatter<thread_probe> : fmt::formatter<std::string_view> {
    auto format(const thread_probe&, fmt::format_context& ctx) const {
        formatted_on_worker.store(std::this_thread::get_id() != main_thread_id, std::memory_order_relaxed);
        return fmt::formatter<std::string_view>::format("probe", ctx);
    }
};

TEST_CASE("Async Print Formats On Worker", "[Utility.AsyncPrint]") {
    suanpan::disable_async_print();

    main_thread_id = std::this_thread::get_id();
    formatted_on_worker.store(false, std::memory_order_relaxed);

    std::ostringstream stream;

    {
        const suanpan::AsyncPrintGuard guard;
        suanpan::print_message(stream, fmt::text_style{}, "{}", thread_probe{});
    }

    REQUIRE(formatted_on_worker.load(std::memory_order_relaxed));
    REQUIRE("probe" == stream.str());
}

TEST_CASE("Async Print Preserves Order", "[Utility.AsyncPrint]") {
    suanpan::disable_async_print();

    std::ostringstream buffer;

    {
        const suanpan::AsyncPrintGuard guard;
        REQUIRE(suanpan::is_async_print_enabled());
        for(int i = 0; i < 1000; ++i) suanpan::print_message(buffer, fmt::text_style{}, "{:04d}\n", i);
    }

    REQUIRE_FALSE(suanpan::is_async_print_enabled());

    std::string expected;
    for(int i = 0; i < 1000; ++i) expected += fmt::format("{:04d}\n", i);
    REQUIRE(expected == buffer.str());
}

TEST_CASE("Async Print Flush", "[Utility.AsyncPrint]") {
    suanpan::disable_async_print();

    std::ostringstream buffer;

    suanpan::print_message(buffer, fmt::text_style{}, "sync\n");
    REQUIRE("sync\n" == buffer.str());

    suanpan::enable_async_print();
    suanpan::print_message(buffer, fmt::text_style{}, "async\n");
    suanpan::flush_print(buffer);
    REQUIRE("sync\nasync\n" == buffer.str());

    suanpan::disable_async_print();
}

TEST_CASE("Async Print Flushes Stream", "[Utility.AsyncPrint]") {
    suanpan::disable_async_print();

    counting_buf buffer;
    std::ostream stream(&buffer);

    suanpan::print_message(stream, fmt::text_style{}, "[====]\r");
    suanpan::flush_print(stream);
    REQUIRE(1 == buffer.sync_count);
    REQUIRE("[====]\r" == buffer.str());

    {
        const suanpan::AsyncPrintGuard guard;
        suanpan::print_message(stream, fmt::text_style{}, "async\r");
        suanpan::flush_print(stream);
        REQUIRE(2 == buffer.sync_count);
    }

    REQUIRE("[====]\rasync\r" == buffer.str());
}

TEST_CASE("Nested Async Print Guards", "[Utility.AsyncPrint]") {
    suanpan::disable_async_print();

    counting_buf buffer;
    std::ostream stream(&buffer);

    {
        const suanpan::AsyncPrintGuard outer;
        {
            const suanpan::AsyncPrintGuard inner;
            suanpan::print_message(stream, fmt::text_style{}, "inner\n");
        }
        REQUIRE(suanpan::is_async_print_enabled());
        suanpan::print_message(stream, fmt::text_style{}, "outer\n");
    }

    REQUIRE_FALSE(suanpan::is_async_print_enabled());
    REQUIRE("inner\nouter\n" == buffer.str());
}

TEST_CASE("Nested Async Print Guard Drains Inner Scope", "[Utility.AsyncPrint]") {
    suanpan::disable_async_print();

    std::ostringstream stream;

    {
        const suanpan::AsyncPrintGuard outer;
        {
            const suanpan::AsyncPrintGuard inner;
            suanpan::print_message(stream, fmt::text_style{}, "inner\n");
        }
        REQUIRE("inner\n" == stream.str());
    }
}

TEST_CASE("Async Print Disable With Concurrent Producers", "[Utility.AsyncPrint]") {
    suanpan::disable_async_print();

    counting_buf buffer;
    std::ostream stream(&buffer);
    constexpr auto producer_count = 4;
    constexpr auto message_count = 100;

    std::atomic go{false};
    std::vector<std::thread> producers;

    {
        const suanpan::AsyncPrintGuard guard;
        for(auto t = 0; t < producer_count; ++t)
            producers.emplace_back([&, t] {
                while(!go.load(std::memory_order_acquire)) {}
                for(auto i = 0; i < message_count; ++i) suanpan::print_message(stream, fmt::text_style{}, "{:d}:{:d}\n", t, i);
            });

        go.store(true, std::memory_order_release);
        // Let the workers submit some messages asynchronously before the guard
        // disables the mode and drains the queue.
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    for(auto& producer : producers) producer.join();

    const auto output = buffer.str();
    for(auto t = 0; t < producer_count; ++t)
        for(auto i = 0; i < message_count; ++i) REQUIRE(std::string::npos != output.find(fmt::format("{:d}:{:d}\n", t, i)));
}

TEST_CASE("Async Print Via Logging Helpers", "[Utility.AsyncPrint]") {
    suanpan::disable_async_print();

    const auto print_backup = SUANPAN_PRINT;
    const auto colour_backup = SUANPAN_COLOR;
    SUANPAN_PRINT = true;
    SUANPAN_COLOR = false;

    auto* buffer_backup = SUANPAN_COUT.rdbuf();
    std::ostringstream buffer;
    SUANPAN_COUT.rdbuf(buffer.rdbuf());

    {
        const suanpan::AsyncPrintGuard guard;
        for(int i = 0; i < 10; ++i) suanpan_info("info {:d}\n", i);
        suanpan_highlight("highlight\n");
    }

    SUANPAN_COUT.rdbuf(buffer_backup);
    SUANPAN_PRINT = print_backup;
    SUANPAN_COLOR = colour_backup;

    const auto output = buffer.str();
    auto pos = std::string::size_type{0};
    for(int i = 0; i < 10; ++i) {
        const auto expected = fmt::format("info {:d}\n", i);
        const auto found = output.find(expected, pos);
        REQUIRE(std::string::npos != found);
        pos = found + expected.size();
    }
    REQUIRE(std::string::npos != output.find("highlight\n", pos));
}

//! @}
