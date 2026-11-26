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
 * @file TestHistory.cpp
 * @brief Tests for loading the CLI command history.
 *
 * @author tlc
 * @date 26/11/2026
 * @version 0.1.0
 * @addtogroup Utility
 * @{
 */

#include "CatchHeader.h"

#include <linenoise/linenoise.hpp>
#include <suanPan.h>

TEST_CASE("Load History Keeps Recent Entries", "[Utility.History]") {
    constexpr auto limit = 64;
    constexpr auto count = 10000;

    const auto path = fs::temp_directory_path() / ("suanpan-history-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    {
        std::ofstream file(path, std::ios::binary);
        REQUIRE(file.is_open());
        for(auto I = 0; I < count; ++I) file << "command " << I << '\n';
        file << "command " << count - 1 << '\n'; // duplicate, must be ignored
    }

    linenoise::detail::history.clear();
    REQUIRE(linenoise::SetHistoryMaxLen(limit));
    REQUIRE(linenoise::LoadHistory(path.generic_string().c_str()));
    fs::remove(path);

    const auto& history = linenoise::GetHistory();
    REQUIRE(history.size() == limit);
    REQUIRE(history.front() == "command " + std::to_string(count - limit));
    REQUIRE(history.back() == "command " + std::to_string(count - 1));
}

TEST_CASE("Load History Appends To Existing Entries", "[Utility.History]") {
    const auto path = fs::temp_directory_path() / ("suanpan-history-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count() + 1));
    {
        std::ofstream file(path, std::ios::binary);
        REQUIRE(file.is_open());
        file << "a\na\nb\nc\nc\nd\n";
    }

    linenoise::detail::history.clear();
    REQUIRE(linenoise::SetHistoryMaxLen(4));
    REQUIRE(linenoise::AddHistory("seed"));
    REQUIRE(linenoise::LoadHistory(path.generic_string().c_str()));
    fs::remove(path);

    const auto& history = linenoise::GetHistory();
    REQUIRE(history.size() == 4);
    REQUIRE(history[0] == "a");
    REQUIRE(history[1] == "b");
    REQUIRE(history[2] == "c");
    REQUIRE(history[3] == "d");
}

//! @}
