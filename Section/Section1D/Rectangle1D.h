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
 * @class Rectangle1D
 * @brief A Rectangle1D class.
 * @author tlc
 * @date 30/10/2017
 * @version 0.1.0
 * @file Rectangle1D.h
 * @addtogroup Section-1D
 * @ingroup Section
 * @{
 */

#ifndef RECTANGLE1D_H
#define RECTANGLE1D_H

#include <Section/Section1D/Section1D.h>

class Rectangle1D : public Section1D {
    const double width, height;

public:
    explicit Rectangle1D(unsigned, // tag
                         double,   // width
                         double,   // height
                         unsigned  // material tag
    );

    unique_ptr<Section> get_copy() override;

    void print() override;
};

#endif

//! @}
