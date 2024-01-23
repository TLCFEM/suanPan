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
 * @class MPF
 * @brief The MPF class.
 *
 * The Menegotto--Pinto steel model.
 *
 * validated @ 11/08/2019 by tlc
 *
 * @author tlc
 * @date 11/08/2019
 * @version 1.0.0
 * @file MPF.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef MPF_H
#define MPF_H

#include <Material/Material1D/Material1D.h>

struct DataMPF {
    const double elastic_modulus;    // elastic modulus
    const double hardening_ratio;    // hardening ratio
    const double yield_stress;       // yield stress
    const double R0, A1, A2, A3, A4; // model parameters

    const bool isotropic_hardening; // isotropic hardening switch
    const bool constant_radius;     // constant radius switch

    const double yield_strain = yield_stress / elastic_modulus; // yield strain
};

class MPF final : protected DataMPF, public Material1D {
public:
    MPF(
        unsigned,      // tag
        double,        // elastic modulus
        double,        // yield stress
        double = .05,  // hardening ratio
        double = 20.,  // R0
        double = 18.5, // A1
        double = .15,  // A2
        double = .01,  // A3
        double = 7.,   // A4
        bool = false,  // isotropic hardening switch
        bool = false,  // constant radius switch
        double = 0.    // density
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
