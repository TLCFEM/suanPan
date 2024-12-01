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
/**
 * @class CDP
 * @brief The CDP class.
 *
 * A 3D concrete material model that supports stiffness degradation.
 *
 * This model uses the original backbones.
 *
 * algorithm verified at 29 April 2019 by tlc
 *
 * References:
 *     1. A Plastic-Damage Model for Concrete.
 *     https://doi.org/10.1016/0020-7683(89)90050-4
 *     2. Plastic-Damage Model for Cyclic Loading of Concrete Structures.
 *     https://doi.org/10.1061/(ASCE)0733-9399(1998)124:8(892)
 *     3. A Plastic-Damage Concrete Model for Earthquake Analysis of Dams.
 *     https://doi.org/10.1002/(SICI)1096-9845(199809)27:9<937::AID-EQE764>3.0.CO;2-5
 *     4. A Return-Mapping Algorithm for Plastic-Damage Models: 3-D and Plane Stress Formulation.
 *     https://doi.org/10.1002/1097-0207(20010120)50:2<487::AID-NME44>3.0.CO;2-N
 *
 * @author tlc
 * @date 29/04/2019
 * @version 1.0.0
 * @file CDP.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef CDP_H
#define CDP_H

#include "NonlinearCDP.h"

class CDP final : public NonlinearCDP {
    const double a_t, cb_t, f_t;
    const double a_c, cb_c, f_c;

    [[nodiscard]] vec6 compute_tension_backbone(double) const override;
    [[nodiscard]] vec6 compute_compression_backbone(double) const override;

public:
    explicit CDP(
        unsigned = 0,     // tag
        double = 3E4,     // elastic modulus
        double = .2,      // poissons ratio
        double = 3.,      // crack stress (+)
        double = 30.,     // crush stress (-)
        double = 1E-3,    // normalized crack energy (+)
        double = 1E-1,    // normalized crush energy (+)
        double = .8,      // hardening after crack stress a_t
        double = 4.,      // hardening after crush stress a_c
        double = .6,      // reference damage factor at half crack stress
        double = .6,      // reference damage factor at crush stress
        double = .2,      // dilatancy parameter
        double = 1.16,    // biaxial compression strength ratio
        double = .5,      // stiffness recovery
        double = 2400E-12 // density
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
