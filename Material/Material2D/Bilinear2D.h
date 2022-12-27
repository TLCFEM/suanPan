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
 * @class Bilinear2D
 * @brief A Bilinear2D class.
 * @author tlc
 * @date 04/07/2018
 * @version 0.2.0
 * @file Bilinear2D.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef BILINEAR2D_H
#define BILINEAR2D_H

#include <Material/Material2D/Material2D.h>
#include <Material/Material3D/vonMises/BilinearJ2.h>

class Bilinear2D final : public Material2D {
    static const uvec F1, F2;

    const double elastic_modulus;
    const double poissons_ratio;

    vec current_full_strain;
    vec trial_full_strain;

    BilinearJ2 base;

    static mat form_stiffness(const mat&);

public:
    Bilinear2D(unsigned,                 // tag
               double,                   // elastic modulus
               double,                   // poisson's ratio
               double,                   // initial yield stress
               double,                   // hardening ratio
               double = 1.,              // isotropic (1.0) / kinematic (0.0) hardening factor
               PlaneType = PlaneType::S, // plane stress or plane strain
               double = 0.               // density
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
