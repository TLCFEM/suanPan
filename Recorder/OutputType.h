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

#ifndef OUTPUTTYPE_H
#define OUTPUTTYPE_H

#include <magic_enum/magic_enum.hpp>

enum class OutputType {
    // stress
    S,
    S11,
    S12,
    S13,
    S22,
    S23,
    S33,

    // strain
    E,
    E11,
    E12,
    E13,
    E22,
    E23,
    E33,

    // elastic strain
    EE,
    EE11,
    EE12,
    EE13,
    EE22,
    EE23,
    EE33,

    // plastic strain
    PE,
    PE11,
    PE12,
    PE13,
    PE22,
    PE23,
    PE33,

    // principal stress
    SP,
    SP1,
    SP2,
    SP3,

    // principal strain
    EP,
    EP1,
    EP2,
    EP3,

    // principal elastic strain
    EEP,
    EEP1,
    EEP2,
    EEP3,

    // principal plastic strain
    PEP,
    PEP1,
    PEP2,
    PEP3,

    // displacement
    U,
    U1,
    U2,
    U3,
    U4,
    U5,
    U6,
    UR1,
    UR2,
    UR3,

    // velocity
    V,
    V1,
    V2,
    V3,
    V4,
    V5,
    V6,
    VR1,
    VR2,
    VR3,

    // acceleration
    A,
    A1,
    A2,
    A3,
    A4,
    A5,
    A6,
    AR1,
    AR2,
    AR3,

    // momentum
    MM,
    MM1,
    MM2,
    MM3,
    MM4,
    MM5,
    MM6,
    MMR1,
    MMR2,
    MMR3,

    // reaction force
    RF,
    RF1,
    RF2,
    RF3,
    RF4,
    RF5,
    RF6,
    RM1,
    RM2,
    RM3,

    // damping force
    DF,
    DF1,
    DF2,
    DF3,
    DF4,
    DF5,
    DF6,
    DM1,
    DM2,
    DM3,

    // inertial force
    IF,
    IF1,
    IF2,
    IF3,
    IF4,
    IF5,
    IF6,
    IM1,
    IM2,
    IM3,

    // global damping force
    GDF,
    GDF1,
    GDF2,
    GDF3,
    GDF4,
    GDF5,
    GDF6,
    GDM1,
    GDM2,
    GDM3,

    // global inertial force
    GIF,
    GIF1,
    GIF2,
    GIF3,
    GIF4,
    GIF5,
    GIF6,
    GIM1,
    GIM2,
    GIM3,

    // hydrostatic pressure
    HYDRO,

    // von Mises stress
    MISES,

    // equivalent strain
    EEQ,
    // equivalent elastic strain
    EEEQ,
    // equivalent plastic strain
    PEEQ,

    // kinetic energy
    KE,
    // strain energy
    SE,
    // viscous energy
    VE,
    // nonviscous energy
    NVE,

    // stiffness
    K,
    // mass
    M,

    // damper strain in Maxwell/Kelvin model
    ED,
    // damper strain rate in Maxwell/Kelvin model
    VD,
    // damper force in Maxwell/Kelvin model
    SD,
    // spring strain in Maxwell/Kelvin model
    ES,
    // spring strain rate in Maxwell/Kelvin model
    VS,
    // spring force in Maxwell/Kelvin model
    SS,

    // damage variable in phase field models
    DAMAGE,
    // tensile damage in material models
    DT,
    // compressive damage in material models
    DC,

    // pore pressure
    PP,

    // volume fraction
    VF,

    // history vector
    HIST,

    // local iteration in Maxwell model
    LITR,

    // yield flag, 1 for yield, 0 for not yield
    YF,

    // beam end deformation
    BEAME,
    // beam end force
    BEAMS,

    // amplitude
    AMP,
    NL
};

template<> struct magic_enum::customize::enum_range<OutputType> {
    static constexpr int min = 0;
    static constexpr int max = 512;
};

constexpr std::string_view to_name(const OutputType L) { return magic_enum::enum_name(L); }

constexpr OutputType to_token(const std::string_view L) { return magic_enum::enum_cast<OutputType>(L).value_or(OutputType::NL); }

std::string to_category(OutputType);

int to_index(OutputType);

#endif
