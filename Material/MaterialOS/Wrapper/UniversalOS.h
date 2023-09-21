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
 * @class UniversalOS
 * @brief A UniversalOS class.
 * @author tlc
 * @date 21/09/2023
 * @version 0.1.0
 * @file UniversalOS.h
 * @addtogroup Material-OS
 * @{
 */

#ifndef UniversalOS_H
#define UniversalOS_H

#include <Material/Material3D/Wrapper/StressWrapper.h>

class UniversalOS : public StressWrapper {
public:
    UniversalOS(unsigned, // tag
                unsigned, // 3D material tag
                unsigned, // max iteration
                uvec&&,
                uvec&&
    );

    void print() override;
};

class OS14 final : public UniversalOS {
public:
    OS14(unsigned, // tag
         unsigned, // 3D material tag
         unsigned  // max iteration
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
