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
 * @class Box2D
 * @brief A Box2D class.
 * @author tlc
 * @date 19/03/2019
 * @version 0.1.0
 * @file Box2D.h
 * @addtogroup Section-2D
 * @ingroup Section
 * @{
 */

#ifndef BOX2D_H
#define BOX2D_H

#include <Section/Section2D/ISection2D.h>

class Box2D final : public ISection2D {
public:
    Box2D(
        unsigned,     // tag
        double,       // width
        double,       // height
        double,       // thickness
        unsigned,     // material tag
        unsigned = 6, // number of integration points
        double = 0.   // eccentricity
    );
    Box2D(
        unsigned,     // tag
        vec&&,        // dimension
        unsigned,     // material tag
        unsigned = 6, // number of integration points
        double = 0.   // eccentricity
    );

    unique_ptr<Section> get_copy() override;

    void print() override;
};

#endif

//! @}
