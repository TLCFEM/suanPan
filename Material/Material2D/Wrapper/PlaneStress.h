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
 * @class PlaneStress
 * @brief A PlaneStress class.
 * @author tlc
 * @date 04/10/2017
 * @version 0.1.0
 * @file PlaneStress.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef PLANESTRESS_H
#define PLANESTRESS_H

#include <Material/Material3D/Wrapper/StressWrapper.h>
#include <Toolbox/ResourceHolder.h>

class PlaneStress final : public StressWrapper {
public:
    PlaneStress(
        unsigned, // tag
        unsigned, // 3D material tag
        unsigned  // max iteration
    );

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
