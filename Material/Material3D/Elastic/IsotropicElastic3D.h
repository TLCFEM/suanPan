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
 * @class IsotropicElastic3D
 * @brief The IsotropicElastic3D class defines a isotropic elastic material for 3-D
 * problems.
 *
 * The Young's modulus is stored in `elastic_modulus`. The Poisson's ratio is
 * stored in `poissons_ratio`.
 *
 * algorithm verified at 24 April 2019 by tlc
 *
 * @author tlc
 * @date 24/04/2019
 * @version 1.0.0
 * @file IsotropicElastic3D.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef ISOTROPICELASTIC3D_H
#define ISOTROPICELASTIC3D_H

#include <Material/Material3D/Material3D.h>

struct DataIsotropicElastic3D {
    double elastic_modulus; // elastic modulus
    double poissons_ratio;  // poissons ratio
};

class IsotropicElastic3D final : public DataIsotropicElastic3D, public Material3D {
public:
    IsotropicElastic3D(unsigned,   // tag
                       double,     // elastic modulus
                       double,     // poissons ratio
                       double = 0. // density
        );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
