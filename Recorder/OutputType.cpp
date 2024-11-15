/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

const char* to_category(const OutputType L) {
    if(OutputType::S == L) return "S";
    if(OutputType::S11 == L) return "S";
    if(OutputType::S12 == L) return "S";
    if(OutputType::S13 == L) return "S";
    if(OutputType::S22 == L) return "S";
    if(OutputType::S23 == L) return "S";
    if(OutputType::S33 == L) return "S";
    if(OutputType::E == L) return "E";
    if(OutputType::E11 == L) return "E";
    if(OutputType::E12 == L) return "E";
    if(OutputType::E13 == L) return "E";
    if(OutputType::E22 == L) return "E";
    if(OutputType::E23 == L) return "E";
    if(OutputType::E33 == L) return "E";
    if(OutputType::EE == L) return "EE";
    if(OutputType::EE11 == L) return "EE";
    if(OutputType::EE12 == L) return "EE";
    if(OutputType::EE13 == L) return "EE";
    if(OutputType::EE22 == L) return "EE";
    if(OutputType::EE23 == L) return "EE";
    if(OutputType::EE33 == L) return "EE";
    if(OutputType::PE == L) return "PE";
    if(OutputType::PE11 == L) return "PE";
    if(OutputType::PE12 == L) return "PE";
    if(OutputType::PE13 == L) return "PE";
    if(OutputType::PE22 == L) return "PE";
    if(OutputType::PE23 == L) return "PE";
    if(OutputType::PE33 == L) return "PE";
    if(OutputType::SP == L) return "SP";
    if(OutputType::SP1 == L) return "SP";
    if(OutputType::SP2 == L) return "SP";
    if(OutputType::SP3 == L) return "SP";
    if(OutputType::EP == L) return "EP";
    if(OutputType::EP1 == L) return "EP";
    if(OutputType::EP2 == L) return "EP";
    if(OutputType::EP3 == L) return "EP";
    if(OutputType::EEP == L) return "EEP";
    if(OutputType::EEP1 == L) return "EEP";
    if(OutputType::EEP2 == L) return "EEP";
    if(OutputType::EEP3 == L) return "EEP";
    if(OutputType::PEP == L) return "PEP";
    if(OutputType::PEP1 == L) return "PEP";
    if(OutputType::PEP2 == L) return "PEP";
    if(OutputType::PEP3 == L) return "PEP";
    if(OutputType::U == L) return "U";
    if(OutputType::U1 == L) return "U";
    if(OutputType::U2 == L) return "U";
    if(OutputType::U3 == L) return "U";
    if(OutputType::U4 == L) return "U";
    if(OutputType::U5 == L) return "U";
    if(OutputType::U6 == L) return "U";
    if(OutputType::UR1 == L) return "U";
    if(OutputType::UR2 == L) return "U";
    if(OutputType::UR3 == L) return "U";
    if(OutputType::V == L) return "V";
    if(OutputType::V1 == L) return "V";
    if(OutputType::V2 == L) return "V";
    if(OutputType::V3 == L) return "V";
    if(OutputType::V4 == L) return "V";
    if(OutputType::V5 == L) return "V";
    if(OutputType::V6 == L) return "V";
    if(OutputType::VR1 == L) return "V";
    if(OutputType::VR2 == L) return "V";
    if(OutputType::VR3 == L) return "V";
    if(OutputType::A == L) return "A";
    if(OutputType::A1 == L) return "A";
    if(OutputType::A2 == L) return "A";
    if(OutputType::A3 == L) return "A";
    if(OutputType::A4 == L) return "A";
    if(OutputType::A5 == L) return "A";
    if(OutputType::A6 == L) return "A";
    if(OutputType::AR1 == L) return "A";
    if(OutputType::AR2 == L) return "A";
    if(OutputType::AR3 == L) return "A";
    if(OutputType::MM == L) return "MM";
    if(OutputType::MM1 == L) return "MM";
    if(OutputType::MM2 == L) return "MM";
    if(OutputType::MM3 == L) return "MM";
    if(OutputType::MM4 == L) return "MM";
    if(OutputType::MM5 == L) return "MM";
    if(OutputType::MM6 == L) return "MM";
    if(OutputType::MMR1 == L) return "MM";
    if(OutputType::MMR2 == L) return "MM";
    if(OutputType::MMR3 == L) return "MM";
    if(OutputType::RF == L) return "RF";
    if(OutputType::RF1 == L) return "RF";
    if(OutputType::RF2 == L) return "RF";
    if(OutputType::RF3 == L) return "RF";
    if(OutputType::RF4 == L) return "RF";
    if(OutputType::RF5 == L) return "RF";
    if(OutputType::RF6 == L) return "RF";
    if(OutputType::RM1 == L) return "RF";
    if(OutputType::RM2 == L) return "RF";
    if(OutputType::RM3 == L) return "RF";
    if(OutputType::DF == L) return "DF";
    if(OutputType::DF1 == L) return "DF";
    if(OutputType::DF2 == L) return "DF";
    if(OutputType::DF3 == L) return "DF";
    if(OutputType::DF4 == L) return "DF";
    if(OutputType::DF5 == L) return "DF";
    if(OutputType::DF6 == L) return "DF";
    if(OutputType::DM1 == L) return "DF";
    if(OutputType::DM2 == L) return "DF";
    if(OutputType::DM3 == L) return "DF";
    if(OutputType::IF == L) return "IF";
    if(OutputType::IF1 == L) return "IF";
    if(OutputType::IF2 == L) return "IF";
    if(OutputType::IF3 == L) return "IF";
    if(OutputType::IF4 == L) return "IF";
    if(OutputType::IF5 == L) return "IF";
    if(OutputType::IF6 == L) return "IF";
    if(OutputType::IM1 == L) return "IF";
    if(OutputType::IM2 == L) return "IF";
    if(OutputType::IM3 == L) return "IF";
    if(OutputType::GDF == L) return "GDF";
    if(OutputType::GDF1 == L) return "GDF";
    if(OutputType::GDF2 == L) return "GDF";
    if(OutputType::GDF3 == L) return "GDF";
    if(OutputType::GDF4 == L) return "GDF";
    if(OutputType::GDF5 == L) return "GDF";
    if(OutputType::GDF6 == L) return "GDF";
    if(OutputType::GDM1 == L) return "GDM";
    if(OutputType::GDM2 == L) return "GDM";
    if(OutputType::GDM3 == L) return "GDM";
    if(OutputType::GIF == L) return "GIF";
    if(OutputType::GIF1 == L) return "GIF";
    if(OutputType::GIF2 == L) return "GIF";
    if(OutputType::GIF3 == L) return "GIF";
    if(OutputType::GIF4 == L) return "GIF";
    if(OutputType::GIF5 == L) return "GIF";
    if(OutputType::GIF6 == L) return "GIF";
    if(OutputType::GIM1 == L) return "GIF";
    if(OutputType::GIM2 == L) return "GIF";
    if(OutputType::GIM3 == L) return "GIF";
    if(OutputType::HYDRO == L) return "HYDRO";
    if(OutputType::MISES == L) return "MISES";
    if(OutputType::EEQ == L) return "EEQ";
    if(OutputType::EEEQ == L) return "EEEQ";
    if(OutputType::PEEQ == L) return "PEEQ";
    if(OutputType::KE == L) return "KE";
    if(OutputType::SE == L) return "SE";
    if(OutputType::VE == L) return "VE";
    if(OutputType::NVE == L) return "NVE";
    if(OutputType::K == L) return "K";
    if(OutputType::M == L) return "M";
    if(OutputType::ED == L) return "ED";
    if(OutputType::VD == L) return "VD";
    if(OutputType::SD == L) return "SD";
    if(OutputType::ES == L) return "ES";
    if(OutputType::VS == L) return "VS";
    if(OutputType::SS == L) return "SS";
    if(OutputType::DAMAGE == L) return "DAMAGE";
    if(OutputType::DT == L) return "DT";
    if(OutputType::DC == L) return "DC";
    if(OutputType::PP == L) return "PP";
    if(OutputType::VF == L) return "VF";
    if(OutputType::HIST == L) return "HIST";
    if(OutputType::LITR == L) return "LITR";
    if(OutputType::YF == L) return "YF";
    if(OutputType::BEAME == L) return "BEAME";
    if(OutputType::BEAMS == L) return "BEAMS";
    if(OutputType::AMP == L) return "AMP";
    if(OutputType::NL == L) return "NL";

    return "NL";
}
