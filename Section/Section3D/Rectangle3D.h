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
 * @class Rectangle3D
 * @brief A Rectangle3D class.
 * @author tlc
 * @date 13/10/2017
 * @version 0.1.0
 * @file Rectangle3D.h
 * @addtogroup Section-3D
 * @ingroup Section
 * @{
 */

#ifndef RECTANGLE3D_H
#define RECTANGLE3D_H

#include <Section/Section3D/Section3D.h>

class Rectangle3D final : public Section3D {
    const double width, height;

    const unsigned int_pt_num;

public:
    Rectangle3D(unsigned,     // tag
                double,       // width
                double,       // height
                unsigned,     // material tag
                unsigned = 6, // number of integration points
                double = 0.,  // eccentricity
                double = 0.   // eccentricity
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    void print() override;
};

#endif

//! @}
