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
/**
 * @class TableGurson
 * @brief A TableGurson material class.
 *
 * @author tlc
 * @date 12/01/2020
 * @version 0.1.0
 * @file TableGurson.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef TABLEGURSON_H
#define TABLEGURSON_H

#include "NonlinearGurson.h"

class TableGurson final : public NonlinearGurson {
    const mat hardening_table;

    [[nodiscard]] vec2 compute_hardening(double) const override;

public:
    TableGurson(
        unsigned,   // tag
        double,     // elastic modulus
        double,     // poisson's ratio
        mat&&,      // table
        double,     // q1
        double,     // q2
        double,     // fn
        double,     // sn
        double,     // en
        double = 0. // density
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
