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
 * @class Box3D
 * @brief A Box3D class.
 * @author tlc
 * @date 19/03/2019
 * @version 0.1.0
 * @file Box3D.h
 * @addtogroup Section-3D
 * @ingroup Section
 * @{
 */

#ifndef BOX3D_H
#define BOX3D_H

#include <Section/Section3D/Section3D.h>

class Box3D final : public Section3D {
    const double width, height, thickness;

    const unsigned int_pt_num;

public:
    Box3D(unsigned,     // tag
          double,       // width
          double,       // height
          double,       // thickness
          unsigned,     // material tag
          unsigned = 6, // number of integration points
          double = 0.,  // eccentricity
          double = 0.   // eccentricity
    );
    Box3D(unsigned,        // tag
          vec&&,           // dimension
          unsigned,        // material tag
          unsigned = 6,    // number of integration points
          vec&& = {0., 0.} // eccentricity
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    void print() override;
};

#endif

//! @}
