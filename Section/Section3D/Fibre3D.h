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
 * @class Fibre3D
 * @brief A Fibre3D class.
 * @author tlc
 * @date 15/09/2023
 * @version 0.1.0
 * @file Fibre3D.h
 * @addtogroup Section-3D
 * @ingroup Section
 * @{
 */

#ifndef FIBRE3D_H
#define FIBRE3D_H

#include <Section/Fibre.h>

class Fibre3D final : public Fibre {
public:
    Fibre3D(unsigned, uvec&&);

    unique_ptr<Section> get_copy() override;
};

#endif

//! @}
