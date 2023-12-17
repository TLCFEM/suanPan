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
 * @class Cell3DOS
 * @brief A Cell3DOS class.
 * @author tlc
 * @date 15/09/2023
 * @version 0.1.0
 * @file Cell3DOS.h
 * @addtogroup Section-OS
 * @ingroup Section
 * @{
 */

#ifndef CELL3DOS_H
#define CELL3DOS_H

#include <Section/SectionOS/SectionOS3D.h>

class Cell3DOS final : public SectionOS3D {
    const double omega, py, pz;

public:
    Cell3DOS(
        unsigned, // tag
        double,   // area
        double,   // sectional coordinate
        double,   // py
        double,   // pz
        unsigned, // material tag
        double,   // eccentricity
        double    // eccentricity
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    void print() override;
};

#endif

//! @}
