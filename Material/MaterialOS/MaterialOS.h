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
 * @class MaterialOS
 * @brief The MaterialOS class.
 *
 * @author tlc
 * @date 15/09/2023
 * @version 0.1.0
 * @file MaterialOS.h
 * @addtogroup Material-OS
 * @ingroup Material
 * @{
 */

#ifndef MATERIALOS_H
#define MATERIALOS_H

#include <Material/Material.h>

using std::vector;

class MaterialOS : public Material {
public:
    MaterialOS(unsigned, // tag
               double    // density
    );
};

#endif

//! @}
