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
 * @class Material1D
 * @brief A Material1D class.
 * @author tlc
 * @date 14/10/2017
 * @version 0.1.0
 * @file Material1D.h
 * @addtogroup Material-1D
 * @ingroup Material
 * @{
 */

#ifndef MATERIAL1D_H
#define MATERIAL1D_H

#include <Material/Material.h>

using std::vector;

class Material1D : public Material {
public:
    Material1D(unsigned, // tag
               double    // density
    );

    const mat& get_trial_stiffness() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
