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
 * @class TSection2D
 * @brief A TSection2D class.
 * @author tlc
 * @date 30/10/2017
 * @version 0.1.0
 * @file TSection2D.h
 * @addtogroup Section-2D
 * @ingroup Section
 * @{
 */

#ifndef TSECTION2D_H
#define TSECTION2D_H

#include <Section/Section2D/ISection2D.h>

class TSection2D final : public ISection2D {
public:
    TSection2D(
        unsigned,     // tag
        double,       // width
        double,       // height
        double,       // width
        double,       // height
        unsigned,     // material tag
        unsigned = 6, // number of integration points
        double = 0.   // eccentricity
    );
    TSection2D(
        unsigned,     // tag
        vec&&,        // dimension
        unsigned,     // material tag
        unsigned = 6, // number of integration points
        double = 0.   // eccentricity
    );

    unique_ptr<Section> get_copy() override;
};

#endif

//! @}
