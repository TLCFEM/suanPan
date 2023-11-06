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
 * @class Cell2D
 * @brief A Cell2D class.
 * @author tlc
 * @date 15/09/2023
 * @version 0.1.1
 * @file Cell2D.h
 * @addtogroup Section-2D
 * @ingroup Section
 * @{
 */

#ifndef CELL2D_H
#define CELL2D_H

#include <Section/Section2D/Section2D.h>

class Cell2D final : public Section2D {
public:
    Cell2D(
        unsigned,   // tag
        double,     // area
        unsigned,   // material tag
        double = 0. // eccentricity
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    void print() override;
};

using Bar2D = Cell2D;

#endif

//! @}
