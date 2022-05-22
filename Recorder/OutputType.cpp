/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
#include <Toolbox/utility.h>

const char* to_char(const OutputType& L) {
    switch(L) {
    case OutputType::SD:
        return "SD";
    case OutputType::ED:
        return "ED";
    case OutputType::VD:
        return "VD";
    case OutputType::SS:
        return "SS";
    case OutputType::ES:
        return "ES";
    case OutputType::VS:
        return "VS";
    case OutputType::S:
        return "S";
    case OutputType::S11:
        return "S11";
    case OutputType::S22:
        return "S22";
    case OutputType::S33:
        return "S33";
    case OutputType::S12:
        return "S12";
    case OutputType::S23:
        return "S23";
    case OutputType::S13:
        return "S13";
    case OutputType::SINT:
        return "SINT";
    case OutputType::HYDRO:
        return "HYDRO";
    case OutputType::PP:
        return "PP";
    case OutputType::E:
        return "E";
    case OutputType::E11:
        return "E11";
    case OutputType::E22:
        return "E22";
    case OutputType::E33:
        return "E33";
    case OutputType::E12:
        return "E12";
    case OutputType::E23:
        return "E23";
    case OutputType::E13:
        return "E13";
    case OutputType::EEQ:
        return "EEQ";
    case OutputType::EINT:
        return "EINT";
    case OutputType::SP:
        return "SP";
    case OutputType::SP1:
        return "SP1";
    case OutputType::SP2:
        return "SP2";
    case OutputType::SP3:
        return "SP3";
    case OutputType::EP:
        return "EP";
    case OutputType::EP1:
        return "EP1";
    case OutputType::EP2:
        return "EP2";
    case OutputType::EP3:
        return "EP3";
    case OutputType::SINV:
        return "SINV";
    case OutputType::MISES:
        return "MISES";
    case OutputType::NMISES:
        return "NMISES";
    case OutputType::TRESC:
        return "TRESC";
    case OutputType::EE:
        return "EE";
    case OutputType::EE11:
        return "EE11";
    case OutputType::EE22:
        return "EE22";
    case OutputType::EE33:
        return "EE33";
    case OutputType::EE12:
        return "EE12";
    case OutputType::EE23:
        return "EE23";
    case OutputType::EE13:
        return "EE13";
    case OutputType::EEP:
        return "EEP";
    case OutputType::EEP1:
        return "EEP1";
    case OutputType::EEP2:
        return "EEP2";
    case OutputType::EEP3:
        return "EEP3";
    case OutputType::EEEQ:
        return "EEEQ";
    case OutputType::PE:
        return "PE";
    case OutputType::PE11:
        return "PE11";
    case OutputType::PE22:
        return "PE22";
    case OutputType::PE33:
        return "PE33";
    case OutputType::PE12:
        return "PE12";
    case OutputType::PE23:
        return "PE23";
    case OutputType::PE13:
        return "PE13";
    case OutputType::PEP:
        return "PEP";
    case OutputType::PEP1:
        return "PEP1";
    case OutputType::PEP2:
        return "PEP2";
    case OutputType::PEP3:
        return "PEP3";
    case OutputType::PEEQ:
        return "PEEQ";

    case OutputType::U:
        return "U";
    case OutputType::UT:
        return "UT";
    case OutputType::UR:
        return "UR";
    case OutputType::U1:
        return "U1";
    case OutputType::U2:
        return "U2";
    case OutputType::U3:
        return "U3";
    case OutputType::UR1:
        return "UR1";
    case OutputType::UR2:
        return "UR2";
    case OutputType::UR3:
        return "UR3";
    case OutputType::U4:
        return "U4";
    case OutputType::U5:
        return "U5";
    case OutputType::U6:
        return "U6";
    case OutputType::V:
        return "V";
    case OutputType::VT:
        return "VT";
    case OutputType::VR:
        return "VR";
    case OutputType::V1:
        return "V1";
    case OutputType::V2:
        return "V2";
    case OutputType::V3:
        return "V3";
    case OutputType::VR1:
        return "VR1";
    case OutputType::VR2:
        return "VR2";
    case OutputType::VR3:
        return "VR3";
    case OutputType::V4:
        return "V4";
    case OutputType::V5:
        return "V5";
    case OutputType::V6:
        return "V6";
    case OutputType::A:
        return "A";
    case OutputType::AT:
        return "AT";
    case OutputType::AR:
        return "AR";
    case OutputType::A1:
        return "A1";
    case OutputType::A2:
        return "A2";
    case OutputType::A3:
        return "A3";
    case OutputType::AR1:
        return "AR1";
    case OutputType::AR2:
        return "AR2";
    case OutputType::AR3:
        return "AR3";
    case OutputType::A4:
        return "A4";
    case OutputType::A5:
        return "A5";
    case OutputType::A6:
        return "A6";

    case OutputType::RF:
        return "RF";
    case OutputType::RT:
        return "RT";
    case OutputType::RF1:
        return "RF1";
    case OutputType::RF2:
        return "RF2";
    case OutputType::RF3:
        return "RF3";
    case OutputType::RF4:
        return "RF4";
    case OutputType::RF5:
        return "RF5";
    case OutputType::RF6:
        return "RF6";
    case OutputType::RM:
        return "RM";
    case OutputType::RM1:
        return "RM1";
    case OutputType::RM2:
        return "RM2";
    case OutputType::RM3:
        return "RM3";
    case OutputType::DF:
        return "DF";
    case OutputType::DF1:
        return "DF1";
    case OutputType::DF2:
        return "DF2";
    case OutputType::DF3:
        return "DF3";
    case OutputType::DF4:
        return "DF4";
    case OutputType::DF5:
        return "DF5";
    case OutputType::DF6:
        return "DF6";
    case OutputType::DM1:
        return "DF4";
    case OutputType::DM2:
        return "DF5";
    case OutputType::DM3:
        return "DF6";
    case OutputType::IF:
        return "IF";
    case OutputType::IF1:
        return "IF1";
    case OutputType::IF2:
        return "IF2";
    case OutputType::IF3:
        return "IF3";
    case OutputType::IF4:
        return "IF4";
    case OutputType::IF5:
        return "IF5";
    case OutputType::IF6:
        return "IF6";
    case OutputType::IM1:
        return "IF4";
    case OutputType::IM2:
        return "IF5";
    case OutputType::IM3:
        return "IF6";
    case OutputType::GDF:
        return "GDF";
    case OutputType::GDF1:
        return "GDF1";
    case OutputType::GDF2:
        return "GDF2";
    case OutputType::GDF3:
        return "GDF3";
    case OutputType::GDF4:
        return "GDF4";
    case OutputType::GDF5:
        return "GDF5";
    case OutputType::GDF6:
        return "GDF6";
    case OutputType::GDM1:
        return "GDF4";
    case OutputType::GDM2:
        return "GDF5";
    case OutputType::GDM3:
        return "GDF6";
    case OutputType::GIF:
        return "GIF";
    case OutputType::GIF1:
        return "GIF1";
    case OutputType::GIF2:
        return "GIF2";
    case OutputType::GIF3:
        return "GIF3";
    case OutputType::GIF4:
        return "GIF4";
    case OutputType::GIF5:
        return "GIF5";
    case OutputType::GIF6:
        return "GIF6";
    case OutputType::GIM1:
        return "GIF4";
    case OutputType::GIM2:
        return "GIF5";
    case OutputType::GIM3:
        return "GIF6";

    case OutputType::DAMAGE:
        return "DAMAGE";
    case OutputType::DT:
        return "DT";
    case OutputType::DC:
        return "DC";
    case OutputType::KAPPAT:
        return "KAPPAT";
    case OutputType::KAPPAC:
        return "KAPPAC";
    case OutputType::KAPPAP:
        return "KAPPAP";
    case OutputType::VF:
        return "VF";

    case OutputType::REBARE:
        return "REBARE";
    case OutputType::REBARS:
        return "REBARS";
    case OutputType::RESULTANT:
        return "RESULTANT";
    case OutputType::AXIAL:
        return "AXIAL";
    case OutputType::SHEAR:
        return "SHEAR";
    case OutputType::MOMENT:
        return "MOMENT";
    case OutputType::TORSION:
        return "TORSION";
    case OutputType::LITR:
        return "LITR";
    case OutputType::K:
        return "STIFFNESS";
    case OutputType::M:
        return "MASS";

    case OutputType::SE:
        return "SE";
    case OutputType::TSE:
        return "TSE";
    case OutputType::CSE:
        return "CSE";
    case OutputType::KE:
        return "KE";
    case OutputType::VE:
        return "VE";
    case OutputType::MM:
        return "MOMENTUM";
    case OutputType::MMX:
        return "MOMENTUMX";
    case OutputType::MMY:
        return "MOMENTUMY";
    case OutputType::MMZ:
        return "MOMENTUMZ";
    case OutputType::MMRX:
        return "MOMENTUMRZ";
    case OutputType::MMRY:
        return "MOMENTUMRY";
    case OutputType::MMRZ:
        return "MOMENTUMRZ";

    case OutputType::HIST:
        return "HIST";

    case OutputType::AMP:
        return "AMPLITUDE";
    case OutputType::YF:
        return "YF";

    case OutputType::NL:
    default:
        return "NL";
    }
}

OutputType to_list(const char* L) {
    if(is_equal(L, "SD")) return OutputType::SD;
    if(is_equal(L, "ED")) return OutputType::ED;
    if(is_equal(L, "VD")) return OutputType::VD;
    if(is_equal(L, "SS")) return OutputType::SS;
    if(is_equal(L, "ES")) return OutputType::ES;
    if(is_equal(L, "VS")) return OutputType::VS;
    if(is_equal(L, "S")) return OutputType::S;
    if(is_equal(L, "S11")) return OutputType::S11;
    if(is_equal(L, "S22")) return OutputType::S22;
    if(is_equal(L, "S33")) return OutputType::S33;
    if(is_equal(L, "S12")) return OutputType::S12;
    if(is_equal(L, "S23")) return OutputType::S23;
    if(is_equal(L, "S13")) return OutputType::S13;
    if(is_equal(L, "SINT")) return OutputType::SINT;
    if(is_equal(L, "HYDRO")) return OutputType::HYDRO;
    if(is_equal(L, "PP")) return OutputType::PP;
    if(is_equal(L, "E")) return OutputType::E;
    if(is_equal(L, "E11")) return OutputType::E11;
    if(is_equal(L, "E22")) return OutputType::E22;
    if(is_equal(L, "E33")) return OutputType::E33;
    if(is_equal(L, "E12")) return OutputType::E12;
    if(is_equal(L, "E23")) return OutputType::E23;
    if(is_equal(L, "E13")) return OutputType::E13;
    if(is_equal(L, "EEQ")) return OutputType::EEQ;
    if(is_equal(L, "EINT")) return OutputType::EINT;
    if(is_equal(L, "SP")) return OutputType::SP;
    if(is_equal(L, "SP1")) return OutputType::SP1;
    if(is_equal(L, "SP2")) return OutputType::SP2;
    if(is_equal(L, "SP3")) return OutputType::SP3;
    if(is_equal(L, "EP")) return OutputType::EP;
    if(is_equal(L, "EP1")) return OutputType::EP1;
    if(is_equal(L, "EP2")) return OutputType::EP2;
    if(is_equal(L, "EP3")) return OutputType::EP3;
    if(is_equal(L, "SINV")) return OutputType::SINV;
    if(is_equal(L, "MISES")) return OutputType::MISES;
    if(is_equal(L, "NMISES")) return OutputType::NMISES;
    if(is_equal(L, "TRESC")) return OutputType::TRESC;

    if(is_equal(L, "EE")) return OutputType::EE;
    if(is_equal(L, "EE11")) return OutputType::EE11;
    if(is_equal(L, "EE22")) return OutputType::EE22;
    if(is_equal(L, "EE33")) return OutputType::EE33;
    if(is_equal(L, "EE12")) return OutputType::EE12;
    if(is_equal(L, "EE23")) return OutputType::EE23;
    if(is_equal(L, "EE13")) return OutputType::EE13;
    if(is_equal(L, "EEP")) return OutputType::EEP;
    if(is_equal(L, "EEP1")) return OutputType::EEP1;
    if(is_equal(L, "EEP2")) return OutputType::EEP2;
    if(is_equal(L, "EEP3")) return OutputType::EEP3;
    if(is_equal(L, "EEEQ")) return OutputType::EEEQ;
    if(is_equal(L, "PE")) return OutputType::PE;
    if(is_equal(L, "PE11")) return OutputType::PE11;
    if(is_equal(L, "PE22")) return OutputType::PE22;
    if(is_equal(L, "PE33")) return OutputType::PE33;
    if(is_equal(L, "PE12")) return OutputType::PE12;
    if(is_equal(L, "PE23")) return OutputType::PE23;
    if(is_equal(L, "PE13")) return OutputType::PE13;
    if(is_equal(L, "PEP")) return OutputType::PEP;
    if(is_equal(L, "PEP1")) return OutputType::PEP1;
    if(is_equal(L, "PEP2")) return OutputType::PEP2;
    if(is_equal(L, "PEP3")) return OutputType::PEP3;
    if(is_equal(L, "PEEQ")) return OutputType::PEEQ;

    if(is_equal(L, "U")) return OutputType::U;
    if(is_equal(L, "UT")) return OutputType::UT;
    if(is_equal(L, "UR")) return OutputType::UR;
    if(is_equal(L, "U1")) return OutputType::U1;
    if(is_equal(L, "U2")) return OutputType::U2;
    if(is_equal(L, "U3")) return OutputType::U3;
    if(is_equal(L, "UR1")) return OutputType::UR1;
    if(is_equal(L, "UR2")) return OutputType::UR2;
    if(is_equal(L, "UR3")) return OutputType::UR3;
    if(is_equal(L, "U4")) return OutputType::U4;
    if(is_equal(L, "U5")) return OutputType::U5;
    if(is_equal(L, "U6")) return OutputType::U6;
    if(is_equal(L, "V")) return OutputType::V;
    if(is_equal(L, "VT")) return OutputType::VT;
    if(is_equal(L, "VR")) return OutputType::VR;
    if(is_equal(L, "V1")) return OutputType::V1;
    if(is_equal(L, "V2")) return OutputType::V2;
    if(is_equal(L, "V3")) return OutputType::V3;
    if(is_equal(L, "VR1")) return OutputType::VR1;
    if(is_equal(L, "VR2")) return OutputType::VR2;
    if(is_equal(L, "VR3")) return OutputType::VR3;
    if(is_equal(L, "V4")) return OutputType::V4;
    if(is_equal(L, "V5")) return OutputType::V5;
    if(is_equal(L, "V6")) return OutputType::V6;
    if(is_equal(L, "A")) return OutputType::A;
    if(is_equal(L, "AT")) return OutputType::AT;
    if(is_equal(L, "AR")) return OutputType::AR;
    if(is_equal(L, "A1")) return OutputType::A1;
    if(is_equal(L, "A2")) return OutputType::A2;
    if(is_equal(L, "A3")) return OutputType::A3;
    if(is_equal(L, "AR1")) return OutputType::AR1;
    if(is_equal(L, "AR2")) return OutputType::AR2;
    if(is_equal(L, "AR3")) return OutputType::AR3;
    if(is_equal(L, "A4")) return OutputType::A4;
    if(is_equal(L, "A5")) return OutputType::A5;
    if(is_equal(L, "A6")) return OutputType::A6;

    if(is_equal(L, "RF")) return OutputType::RF;
    if(is_equal(L, "RT")) return OutputType::RT;
    if(is_equal(L, "RF1")) return OutputType::RF1;
    if(is_equal(L, "RF2")) return OutputType::RF2;
    if(is_equal(L, "RF3")) return OutputType::RF3;
    if(is_equal(L, "RF4")) return OutputType::RF4;
    if(is_equal(L, "RF5")) return OutputType::RF5;
    if(is_equal(L, "RF6")) return OutputType::RF6;
    if(is_equal(L, "RM")) return OutputType::RM;
    if(is_equal(L, "RM1")) return OutputType::RM1;
    if(is_equal(L, "RM2")) return OutputType::RM2;
    if(is_equal(L, "RM3")) return OutputType::RM3;
    if(is_equal(L, "DF")) return OutputType::DF;
    if(is_equal(L, "DF1")) return OutputType::DF1;
    if(is_equal(L, "DF2")) return OutputType::DF2;
    if(is_equal(L, "DF3")) return OutputType::DF3;
    if(is_equal(L, "DF4")) return OutputType::DF4;
    if(is_equal(L, "DF5")) return OutputType::DF5;
    if(is_equal(L, "DF6")) return OutputType::DF6;
    if(is_equal(L, "DM1")) return OutputType::DF4;
    if(is_equal(L, "DM2")) return OutputType::DF5;
    if(is_equal(L, "DM3")) return OutputType::DF6;
    if(is_equal(L, "IF")) return OutputType::IF;
    if(is_equal(L, "IF1")) return OutputType::IF1;
    if(is_equal(L, "IF2")) return OutputType::IF2;
    if(is_equal(L, "IF3")) return OutputType::IF3;
    if(is_equal(L, "IF4")) return OutputType::IF4;
    if(is_equal(L, "IF5")) return OutputType::IF5;
    if(is_equal(L, "IF6")) return OutputType::IF6;
    if(is_equal(L, "IM1")) return OutputType::IF4;
    if(is_equal(L, "IM2")) return OutputType::IF5;
    if(is_equal(L, "IM3")) return OutputType::IF6;
    if(is_equal(L, "GDF")) return OutputType::GDF;
    if(is_equal(L, "GDF1")) return OutputType::GDF1;
    if(is_equal(L, "GDF2")) return OutputType::GDF2;
    if(is_equal(L, "GDF3")) return OutputType::GDF3;
    if(is_equal(L, "GDF4")) return OutputType::GDF4;
    if(is_equal(L, "GDF5")) return OutputType::GDF5;
    if(is_equal(L, "GDF6")) return OutputType::GDF6;
    if(is_equal(L, "GDM1")) return OutputType::GDF4;
    if(is_equal(L, "GDM2")) return OutputType::GDF5;
    if(is_equal(L, "GDM3")) return OutputType::GDF6;
    if(is_equal(L, "GIF")) return OutputType::GIF;
    if(is_equal(L, "GIF1")) return OutputType::GIF1;
    if(is_equal(L, "GIF2")) return OutputType::GIF2;
    if(is_equal(L, "GIF3")) return OutputType::GIF3;
    if(is_equal(L, "GIF4")) return OutputType::GIF4;
    if(is_equal(L, "GIF5")) return OutputType::GIF5;
    if(is_equal(L, "GIF6")) return OutputType::GIF6;
    if(is_equal(L, "GIM1")) return OutputType::GIF4;
    if(is_equal(L, "GIM2")) return OutputType::GIF5;
    if(is_equal(L, "GIM3")) return OutputType::GIF6;

    if(is_equal(L, "DAMAGE")) return OutputType::DAMAGE;
    if(is_equal(L, "DT")) return OutputType::DT;
    if(is_equal(L, "DC")) return OutputType::DC;
    if(is_equal(L, "KAPPAT")) return OutputType::KAPPAT;
    if(is_equal(L, "KAPPAC")) return OutputType::KAPPAC;
    if(is_equal(L, "KAPPAP")) return OutputType::KAPPAP;
    if(is_equal(L, "VF")) return OutputType::VF;

    if(is_equal(L, "REBARE")) return OutputType::REBARE;
    if(is_equal(L, "REBARS")) return OutputType::REBARS;
    if(is_equal(L, "RESULTANT")) return OutputType::RESULTANT;
    if(is_equal(L, "AXIAL")) return OutputType::AXIAL;
    if(is_equal(L, "SHEAR")) return OutputType::SHEAR;
    if(is_equal(L, "MOMENT")) return OutputType::MOMENT;
    if(is_equal(L, "TORSION")) return OutputType::TORSION;

    if(is_equal(L, "LITR")) return OutputType::LITR;
    if(is_equal(L, "STIFFNESS")) return OutputType::K;
    if(is_equal(L, "MASS")) return OutputType::M;
    if(is_equal(L, "K")) return OutputType::K;
    if(is_equal(L, "M")) return OutputType::M;

    if(is_equal(L, "SE")) return OutputType::SE;
    if(is_equal(L, "TSE")) return OutputType::TSE;
    if(is_equal(L, "CSE")) return OutputType::CSE;
    if(is_equal(L, "KE")) return OutputType::KE;
    if(is_equal(L, "VE")) return OutputType::VE;
    if(is_equal(L, "MOMENTUM")) return OutputType::MM;
    if(is_equal(L, "MOMENTUMX")) return OutputType::MMX;
    if(is_equal(L, "MOMENTUMY")) return OutputType::MMY;
    if(is_equal(L, "MOMENTUMZ")) return OutputType::MMZ;
    if(is_equal(L, "MOMENTUMRX")) return OutputType::MMRX;
    if(is_equal(L, "MOMENTUMRY")) return OutputType::MMRY;
    if(is_equal(L, "MOMENTUMRZ")) return OutputType::MMRZ;

    if(is_equal(L, "HIST")) return OutputType::HIST;
    if(is_equal(L, "AMP")) return OutputType::AMP;
    if(is_equal(L, "YF")) return OutputType::YF;

    if(is_equal(L, "NL")) return OutputType::NL;

    return OutputType::NL;
}
