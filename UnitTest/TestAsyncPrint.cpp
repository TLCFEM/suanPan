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
#include <string>
#include <suanPan.h>

TEST_CASE("Async Print Preserves Order", "[Utility.AsyncPrint]") {
    suanpan::disable_async_print();

    std::ostringstream buffer;

    {
        const suanpan::AsyncPrintGuard guard;
        REQUIRE(suanpan::is_async_print_enabled());
        for(int i = 0; i < 1000; ++i) suanpan::print_message(fmt::format("{:04d}\n", i), buffer);
    }

    REQUIRE_FALSE(suanpan::is_async_print_enabled());

    std::string expected;
    for(int i = 0; i < 1000; ++i) expected += fmt::format("{:04d}\n", i);
    REQUIRE(expected == buffer.str());
}

TEST_CASE("Async Print Flush", "[Utility.AsyncPrint]") {
    suanpan::disable_async_print();

    std::ostringstream buffer;

    suanpan::print_message("sync\n", buffer);
    REQUIRE("sync\n" == buffer.str());

    suanpan::enable_async_print();
    suanpan::print_message("async\n", buffer);
    suanpan::flush_print();
    REQUIRE("sync\nasync\n" == buffer.str());

    suanpan::disable_async_print();
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
