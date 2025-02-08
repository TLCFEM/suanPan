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

#include "OutputType.h"

#include <regex>

int to_index(const OutputType config) {
    if(config == OutputType::S11) return 0;
    if(config == OutputType::S22) return 1;
    if(config == OutputType::S33) return 2;
    if(config == OutputType::S12) return 3;
    if(config == OutputType::S23) return 4;
    if(config == OutputType::S13) return 5;

    if(config == OutputType::E11) return 0;
    if(config == OutputType::E22) return 1;
    if(config == OutputType::E33) return 2;
    if(config == OutputType::E12) return 3;
    if(config == OutputType::E23) return 4;
    if(config == OutputType::E13) return 5;

    if(config == OutputType::EE11) return 0;
    if(config == OutputType::EE22) return 1;
    if(config == OutputType::EE33) return 2;
    if(config == OutputType::EE12) return 3;
    if(config == OutputType::EE23) return 4;
    if(config == OutputType::EE13) return 5;

    if(config == OutputType::PE11) return 0;
    if(config == OutputType::PE22) return 1;
    if(config == OutputType::PE33) return 2;
    if(config == OutputType::PE12) return 3;
    if(config == OutputType::PE23) return 4;
    if(config == OutputType::PE13) return 5;

    if(config == OutputType::SP1) return 0;
    if(config == OutputType::SP2) return 1;
    if(config == OutputType::SP3) return 2;

    if(config == OutputType::EP1) return 0;
    if(config == OutputType::EP2) return 1;
    if(config == OutputType::EP3) return 2;

    if(config == OutputType::EEP1) return 0;
    if(config == OutputType::EEP2) return 1;
    if(config == OutputType::EEP3) return 2;

    if(config == OutputType::PEP1) return 0;
    if(config == OutputType::PEP2) return 1;
    if(config == OutputType::PEP3) return 2;

    if(config == OutputType::U1) return 0;
    if(config == OutputType::U2) return 1;
    if(config == OutputType::U3) return 2;
    if(config == OutputType::U4) return 3;
    if(config == OutputType::U5) return 4;
    if(config == OutputType::U6) return 5;
    if(config == OutputType::UR1) return 3;
    if(config == OutputType::UR2) return 4;
    if(config == OutputType::UR3) return 5;

    if(config == OutputType::V1) return 0;
    if(config == OutputType::V2) return 1;
    if(config == OutputType::V3) return 2;
    if(config == OutputType::V4) return 3;
    if(config == OutputType::V5) return 4;
    if(config == OutputType::V6) return 5;
    if(config == OutputType::VR1) return 3;
    if(config == OutputType::VR2) return 4;
    if(config == OutputType::VR3) return 5;

    if(config == OutputType::A1) return 0;
    if(config == OutputType::A2) return 1;
    if(config == OutputType::A3) return 2;
    if(config == OutputType::A4) return 3;
    if(config == OutputType::A5) return 4;
    if(config == OutputType::A6) return 5;
    if(config == OutputType::AR1) return 3;
    if(config == OutputType::AR2) return 4;
    if(config == OutputType::AR3) return 5;

    if(config == OutputType::RF1) return 0;
    if(config == OutputType::RF2) return 1;
    if(config == OutputType::RF3) return 2;
    if(config == OutputType::RF4) return 3;
    if(config == OutputType::RF5) return 4;
    if(config == OutputType::RF6) return 5;
    if(config == OutputType::RM1) return 3;
    if(config == OutputType::RM2) return 4;
    if(config == OutputType::RM3) return 5;

    if(config == OutputType::DF1) return 0;
    if(config == OutputType::DF2) return 1;
    if(config == OutputType::DF3) return 2;
    if(config == OutputType::DF4) return 3;
    if(config == OutputType::DF5) return 4;
    if(config == OutputType::DF6) return 5;
    if(config == OutputType::DM1) return 3;
    if(config == OutputType::DM2) return 4;
    if(config == OutputType::DM3) return 5;

    if(config == OutputType::IF1) return 0;
    if(config == OutputType::IF2) return 1;
    if(config == OutputType::IF3) return 2;
    if(config == OutputType::IF4) return 3;
    if(config == OutputType::IF5) return 4;
    if(config == OutputType::IF6) return 5;
    if(config == OutputType::IM1) return 3;
    if(config == OutputType::IM2) return 4;
    if(config == OutputType::IM3) return 5;

    return 0;
}

std::string to_category(const OutputType L) {
    auto result = std::regex_replace(std::string(to_name(L)), std::regex(R"(\d)"), "");
    if(result == "UR") return "U";
    if(result == "VR") return "V";
    if(result == "AR") return "A";
    if(result == "MMR") return "MM";
    if(result == "RM") return "RF";
    if(result == "DM") return "DF";
    if(result == "IM") return "IF";
    if(result == "GDM") return "GDF";
    if(result == "GIM") return "GIF";
    return result;
}
