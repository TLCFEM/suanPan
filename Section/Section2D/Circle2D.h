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
 * @class Circle2D
 * @brief A Circle2D class.
 * @author tlc
 * @date 15/09/2017
 * @version 0.1.0
 * @file Circle2D.h
 * @addtogroup Section-2D
 * @ingroup Section
 * @{
 */

#ifndef CIRCLE2D_H
#define CIRCLE2D_H

#include <Section/Section2D/Section2D.h>

class Circle2D final : public Section2D {
    const double radius;

    const unsigned int_pt_num;

public:
    Circle2D(
        unsigned,     // tag
        double,       // radius
        unsigned,     // material tag
        unsigned = 6, // number of integration points
        double = 0.   // eccentricity
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    void print() override;
};

#endif

//! @}
