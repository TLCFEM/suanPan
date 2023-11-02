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
 * @class BilinearElastic1D
 * @brief A BilinearElastic1D material class.
 * @author tlc
 * @date 08/08/2017
 * @version 0.1.0
 * @file BilinearElastic1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef BILINEARELASTIC1D_H
#define BILINEARELASTIC1D_H

#include <Material/Material1D/Material1D.h>

struct DataBilinearElastic1D {
    const double elastic_modulus;   // elastic modulus
    const double yield_stress;      // initial yield stress
    const double hardening_modulus; // hardening modulus
    const double radius;            // radius
};

class BilinearElastic1D final : protected DataBilinearElastic1D, public Material1D {
    const double yield_strain = yield_stress / elastic_modulus;

public:
    BilinearElastic1D(
        unsigned, // tag
        double,   // elastic modulus
        double,   // initial yield stress
        double,   // hardening ratio
        double,   // radius
        double    // density
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
