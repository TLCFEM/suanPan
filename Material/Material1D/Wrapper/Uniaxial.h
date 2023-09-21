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
 * @class Uniaxial
 * @brief A Uniaxial class.
 * @author tlc
 * @date 22/01/2019
 * @version 0.1.0
 * @file Uniaxial.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef UNIAXIAL_H
#define UNIAXIAL_H

#include <Material/Material3D/Wrapper/StressWrapper.h>

class Uniaxial final : public StressWrapper {
public:
    Uniaxial(unsigned, // tag
             unsigned, // 3D material tag
             unsigned  // max iteration
    );

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
