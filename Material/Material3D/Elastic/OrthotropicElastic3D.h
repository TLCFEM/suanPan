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
 * @class OrthotropicElastic3D
 * @brief The OrthotropicElastic3D class defines an orthotropic elastic material for 3-D
 * problems.
 *
 * algorithm verified at 24 April 2019 by tlc
 *
 * @author tlc
 * @date 24/04/2019
 * @version 1.0.0
 * @file OrthotropicElastic3D.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef ORTHOTROPICELASTIC3D_H
#define ORTHOTROPICELASTIC3D_H

#include <Material/Material3D/Material3D.h>

class OrthotropicElastic3D final : public Material3D {
    const vec modulus, poissons_ratio;

public:
    OrthotropicElastic3D(unsigned,   // tag
                         vec&&,      // elastic modulus
                         vec&&,      // poissons ratio
                         double = 0. // density
        );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
