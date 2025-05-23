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
 * @class TrussSection
 * @brief A TrussSection class.
 * @author tlc
 * @date 30/10/2017
 * @version 0.1.0
 * @file TrussSection.h
 * @addtogroup Section-1D
 * @ingroup Section
 * @{
 */

#ifndef TRUSSSECTION_H
#define TRUSSSECTION_H

#include <Section/Section1D/Section1D.h>

class TrussSection final : public Section1D {
public:
    explicit TrussSection(
        unsigned, // tag
        double,   // area
        unsigned  // material tag
    );

    unique_ptr<Section> get_copy() override;

    void print() override;
};

#endif

//! @}
