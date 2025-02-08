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
 * @class HSection2D
 * @brief A HSection2D class.
 * @author tlc
 * @date 13/10/2017
 * @version 0.1.0
 * @file HSection2D.h
 * @addtogroup Section-2D
 * @ingroup Section
 * @{
 */

#ifndef HSECTION2D_H
#define HSECTION2D_H

#include <Section/Section2D/Section2D.h>

class HSection2D final : public Section2D {
    const double left_flange_height, left_flange_thickness;
    const double right_flange_height, right_flange_thickness;
    const double web_width, web_thickness;

    const unsigned int_pt_num;

public:
    HSection2D(
        unsigned,     // tag
        double,       // width
        double,       // height
        double,       // width
        double,       // height
        double,       // width
        double,       // height
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
