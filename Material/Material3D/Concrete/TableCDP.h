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
 * @class TableCDP
 * @brief The TableCDP class.
 *
 * A 3D concrete material model that supports stiffness degradation.
 *
 * @author tlc
 * @date 01/02/2020
 * @version 1.0.0
 * @file TableCDP.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef TABLECDP_H
#define TABLECDP_H

#include "NonlinearCDP.h"

class TableCDP final : public NonlinearCDP {
    mat t_table, c_table, dt_table, dc_table;

    [[nodiscard]] vec6 compute_tension_backbone(double) const override;
    [[nodiscard]] vec6 compute_compression_backbone(double) const override;

public:
    TableCDP(
        unsigned,         // tag
        double,           // elastic modulus
        double,           // poissons ratio
        mat&&,            // tension table
        mat&&,            // compression table
        mat&&,            // tension damage table
        mat&&,            // compression damage table
        double = .2,      // dilatancy parameter
        double = 1.16,    // biaxial compression strength ratio
        double = .5,      // stiffness recovery
        double = 2400E-12 // density
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
