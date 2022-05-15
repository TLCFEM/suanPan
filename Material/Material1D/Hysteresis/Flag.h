/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
 * @class Flag
 * @brief A Flag material class.
 * @author tlc
 * @date 14/07/2018
 * @version 0.1.0
 * @file Flag.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef FLAG_H
#define FLAG_H

#include <Material/Material1D/Material1D.h>

struct DataFlag {
    const double elastic_modulus; // elastic modulus

    const double t_hardening_ratio;
    const double t_yield_stress;
    const double t_residual_stress;

    const double c_hardening_ratio;
    const double c_yield_stress;
    const double c_residual_stress;

    const double t_yield_strain = t_yield_stress / elastic_modulus;
    const double t_residual_strain = t_residual_stress / elastic_modulus;
    const double c_yield_strain = c_yield_stress / elastic_modulus;
    const double c_residual_strain = c_residual_stress / elastic_modulus;
};

class Flag final : DataFlag, public Material1D {
    enum class Status {
        NONE,
        TLOAD,
        TUNLOAD,
        TLOW,
        CLOAD,
        CUNLOAD,
        CLOW
    };

    Status trial_status = Status::NONE, current_status = Status::NONE;

public:
    Flag(unsigned, // tag
         double,   // elastic modulus
         double,   // tension initial yield stress
         double,   // tension residual stress
         double,   // tension hardening ratio
         double,   // compression initial yield stress
         double,   // compression residual stress
         double,   // compression hardening ratio
         double    // density
    );
    Flag(unsigned, // tag
         double,   // elastic modulus
         double,   // initial yield stress
         double,   // residual stress
         double,   // hardening ratio
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
