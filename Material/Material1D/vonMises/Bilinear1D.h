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
 * @class Bilinear1D
 * @brief A Bilinear1D material class.
 * @author tlc
 * @date 11/01/2019
 * @version 0.2.0
 * @file Bilinear1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef BILINEAR1D_H
#define BILINEAR1D_H

#include <Material/Material1D/Material1D.h>

struct DataBilinear1D {
    const double elastic_modulus; // elastic modulus
    const double yield_stress;    // initial yield stress
    const double beta;            // isotropic (1.0) / kinematic (0.0) hardening factor
    const double isotropic_modulus;
    const double kinematic_modulus;
};

class Bilinear1D final : protected DataBilinear1D, public Material1D {
public:
    Bilinear1D(
        unsigned,     // tag
        double,       // elastic modulus
        double,       // initial yield stress
        double = .05, // hardening ratio
        double = 1.,  // isotropic/kinematic hardening factor
        double = 0.   // density
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
