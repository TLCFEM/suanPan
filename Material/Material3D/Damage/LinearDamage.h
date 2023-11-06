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
 * @class LinearDamage
 * @brief A LinearDamage material class.
 * @author tlc
 * @date 09/06/2020
 * @version 0.1.0
 * @file LinearDamage.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef LINEARDAMAGE_H
#define LINEARDAMAGE_H

#include "IsotropicDamage.h"

class LinearDamage final : public IsotropicDamage {
    static const double root_two_third;

    const double s_strain, e_strain, e_damage;

    const double slope = (1. - e_damage) / (s_strain - e_strain);

protected:
    void compute_damage() override;

public:
    LinearDamage(
        unsigned, // tag
        unsigned, // mat tag
        double,   // start strain
        double,   // end strain
        double    // end damage value
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
