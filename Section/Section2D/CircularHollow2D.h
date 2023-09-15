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
 * @class CircularHollow2D
 * @brief A CircularHollow2D class.
 * @author tlc
 * @date 14/03/2019
 * @version 0.1.0
 * @file CircularHollow2D.h
 * @addtogroup Section-2D
 * @ingroup Section
 * @{
 */

#ifndef CIRCULARHOLLOW2D_H
#define CIRCULARHOLLOW2D_H

#include <Section/Section2D/Section2D.h>

class CircularHollow2D final : public Section2D {
    const double radius, thickness;

    const unsigned int_pt_num;

public:
    CircularHollow2D(unsigned,      // tag
                     double,        // radius
                     double,        // thickness
                     unsigned,      // material tag
                     unsigned = 10, // number of integration points
                     double = 0.    // eccentricity
    );
    CircularHollow2D(unsigned,      // tag
                     vec&&,         // dimension
                     unsigned,      // material tag
                     unsigned = 10, // number of integration points
                     double = 0.    // eccentricity
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    void print() override;
};

#endif

//! @}
