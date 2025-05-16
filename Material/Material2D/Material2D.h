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
 * @class Material2D
 * @brief The Material2D class defines a isotropic elastic material for plane stress and plane strain problems.
 *
 * The Young's modulus is stored in `elastic_modulus`. The Poisson's ratio is
 * stored in `poissons_ratio`. The `plane_type` labels if it is plane stress or
 * plane strain. The default value `PlaneType::S` represents plane stress.
 * Initializing the object with a `PlaneType::E` value gives a plane strain type
 * response.
 *
 * @author tlc
 * @date 04/10/2017
 * @version 0.1.2
 * @file Material2D.h
 * @addtogroup Material-2D
 * @ingroup Material
 * @{
 */

#ifndef MATERIAL2D_H
#define MATERIAL2D_H

#include <Material/Material.h>

class Material2D : public Material {
public:
    Material2D(
        unsigned,  // tag
        PlaneType, // plane type
        double     // density
    );

    std::vector<vec> record(OutputType) override;
};

#endif

//! @}
