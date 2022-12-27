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
 * @class TSection3D
 * @brief A TSection3D class.
 * @author tlc
 * @date 13/03/2019
 * @version 0.1.0
 * @file TSection3D.h
 * @addtogroup Section-3D
 * @ingroup Section
 * @{
 */

#ifndef TSECTION3D_H
#define TSECTION3D_H

#include <Section/Section3D/ISection3D.h>

class TSection3D final : public ISection3D {
public:
    TSection3D(unsigned,        // tag
               double,          // width
               double,          // height
               double,          // width
               double,          // height
               unsigned,        // material tag
               unsigned = 6,    // number of integration points
               vec&& = {0., 0.} // eccentricity
    );
    TSection3D(unsigned,        // tag
               vec&&,           // dimension
               unsigned,        // material tag
               unsigned = 6,    // number of integration points
               vec&& = {0., 0.} // eccentricity
    );

    unique_ptr<Section> get_copy() override;
};

#endif

//! @}
