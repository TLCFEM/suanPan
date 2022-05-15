/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
 * @class Material3D
 * @brief The Material3D class.
 *
 * @author tlc
 * @date 14/10/2017
 * @version 0.1.0
 * @file Material3D.h
 * @addtogroup Material-3D
 * @ingroup Material
 * @{
 */

#ifndef MATERIAL3D_H
#define MATERIAL3D_H

#include <Material/Material.h>

class Material3D : public Material {
public:
    Material3D(unsigned, // tag
               double    // density
    );

    vector<vec> record(OutputType) override;
};

#endif

//! @}
