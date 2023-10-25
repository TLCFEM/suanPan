/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include <string>

enum class OutputType {
    // damper stress in Maxwell model
    SD,
    // damper strain in Maxwell model
    ED,
    // damper strain rate in Maxwell model
    VD,
    // spring stress in Maxwell model
    SS,
    // spring strain in Maxwell model
    ES,
    // spring strain rate in Maxwell model
    VS,
    // stress
    S,
    S11,
    S22,
    S33,
    S12,
    S23,
    S13,
    // interpolation parameters of the stress field
    SINT,
    HYDRO,
    PP,
    // pore pressure
    // strain
    E,
    E11,
    E22,
    E33,
    E12,
    E23,
    E13,
    // eqv. strain
    EEQ,
    // interpolation parameters of the strain field
    EINT,
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
    // stress invariant
    SINV,
    // Mises stress
    MISES,
    // Mises stress
    NMISES,
    // Tresca stress
    TRESC,
    // elastic strain
    EE,
    EE11,
    EE22,
    EE33,
    EE12,
    EE23,
    EE13,
    // principal elastic strain
    EEP,
    EEP1,
    EEP2,
    EEP3,
    // equivalent elastic strain
    EEEQ,
    // plastic strain
    PE,
    PE11,
    PE22,
    PE33,
    PE12,
    PE23,
    PE13,
    // principal plastic strain
    PEP,
    PEP1,
    PEP2,
    PEP3,
    // eqv. plastic strain
    PEEQ,

    // displacement
    U,
    // translation displacement
    UT,
    // rotation displacement
    UR,
    U1,
    U2,
    U3,
    UR1,
    UR2,
    UR3,
    U4,
    U5,
    U6,
    // velocity
    V,
    VT,
    VR,
    V1,
    V2,
    V3,
    VR1,
    VR2,
    VR3,
    V4,
    V5,
    V6,
    // acceleration
    A,
    AT,
    AR,
    A1,
    A2,
    A3,
    AR1,
    AR2,
    AR3,
    A4,
    A5,
    A6,

    // reaction force
    RF,
    RT,
    RF1,
    RF2,
    RF3,
    RF4,
    RF5,
    RF6,
    // reaction moment
    RM,
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
    // damping force
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
    // inertial force
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

    // damage variable
    DAMAGE,
    DT,
    DC,
    KAPPAT,
    KAPPAC,
    KAPPAP,
    VF,

    RESULTANT,
    // resultant force
    AXIAL,
    // axial force resultant force
    SHEAR,
    // shear force resultant force
    MOMENT,
    // moment resultant force
    TORSION,
    // torsion resultant force

    REBARE,
    REBARS,

    LITR,
    // local iteration counter
    K,
    M,

    // energy
    SE,
    TSE,
    CSE,
    KE,
    VE,
    NVE,
    // momentum
    MM,
    MMX,
    MMY,
    MMZ,
    MMRX,
    MMRY,
    MMRZ,

    // history variable
    HIST,

    // amplitude
    AMP,

    // yield flag
    YF,

    // beam element end forces
    BEAME,
    BEAMS,

    NL
};

const char* to_char(const OutputType&);
OutputType to_list(const char*);
OutputType to_list(const std::string&);

#endif
