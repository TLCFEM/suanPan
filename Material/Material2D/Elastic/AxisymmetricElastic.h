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
 * @class AxisymmetricElastic
 * @brief The AxisymmetricElastic class defines a isotropic elastic material.
 *
 * @author tlc
 * @date 16/02/2019
 * @version 0.1.0
 * @file AxisymmetricElastic.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef AXISTMMETRICELASTIC_H
#define AXISTMMETRICELASTIC_H

#include <Material/Material2D/Material2D.h>

class AxisymmetricElastic final : public Material2D {
    const double elastic_modulus; // elastic modulus
    const double poissons_ratio;  // poissons ratio
public:
    AxisymmetricElastic(unsigned,  // tag
                        double,    // elastic modulus
                        double,    // poissons ratio
                        double = 0 // density
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
