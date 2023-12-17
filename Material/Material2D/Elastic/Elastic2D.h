/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class Elastic2D
 * @brief The Elastic2D class defines a isotropic elastic material for
 * plane stress and
 * plane strain problems.
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
 * @file Elastic2D.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef ELASTIC2D_H
#define ELASTIC2D_H

#include <Material/Material2D/Material2D.h>

class Elastic2D final : public Material2D {
    const double elastic_modulus; // elastic modulus
    const double poissons_ratio;  // poissons ratio
public:
    Elastic2D(
        unsigned,   // tag
        double,     // elastic modulus
        double,     // poissons ratio
        double = 0, // density
        PlaneType = PlaneType::S
    );

    int initialize(const shared_ptr<DomainBase>&) override;
    void initialize_couple(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;

    vector<vec> record(OutputType) override;
};

#endif

//! @}
