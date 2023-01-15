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
/**
 * @class CustomCDP
 * @brief The CustomCDP class.
 *
 * A 3D concrete material model that supports stiffness degradation.
 *
 * This model uses custom backbones.
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
 * @date 16/01/2023
 * @version 1.0.0
 * @file CustomCDP.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef CUSTOMCDP_H
#define CUSTOMCDP_H

#include "NonlinearCDP.h"
#include <Toolbox/Expression.h>

class CustomCDP final : public NonlinearCDP {
    const unsigned t_tag, c_tag;

    shared_ptr<Expression> t_expression, c_expression;

    [[nodiscard]] podarray<double> compute_tension_backbone(double) const override;
    [[nodiscard]] podarray<double> compute_compression_backbone(double) const override;

public:
    CustomCDP(unsigned, // tag
              unsigned, // tension expression tag
              unsigned, // compression expression tag
              double,   // elastic modulus
              double,   // poissons ratio
              double,   // normalized crack energy (+)
              double,   // normalized crush energy (+)
              double,   // dilatancy parameter
              double,   // biaxial compression strength ratio
              double,   // stiffness recovery
              double    // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
