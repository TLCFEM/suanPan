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

#ifndef SOLVERSETTING_HPP
#define SOLVERSETTING_HPP

#include <Toolbox/utility.h>

enum class Precision : std::uint8_t {
    MIXED,
    FULL
};

template<sp_d data_t> struct SolverSetting {
    std::string option{};
    data_t tolerance = std::is_same_v<data_t, float> ? 1E-7f : 1E-14;
    std::uint8_t iterative_refinement = 5;
    Precision precision = Precision::FULL;

    auto set_option(istringstream& command) { option = get_remaining(command); }
    auto set_lis_option(istringstream& command) {
        static constexpr auto max_length = 1024;

        const auto sub_command = get_remaining(command);

        if(sub_command.empty()) {
            option = "-i fgmres -p ilu";
            return;
        }
        if(std::any_of(sub_command.begin(), sub_command.end(), [](const char c) { return !std::isspace(c); })) option = sub_command;
        if(option.length() < max_length) return;

        const auto pos = option.find_last_of(" \t\n\r", max_length);
        option = option.substr(0, pos == std::string::npos ? max_length : pos);
    }
};

#endif
