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
 * @class AFC
 * @brief A AFC material class.
 * @author tlc
 * @date 30/03/2019
 * @version 0.1.0
 * @file AFC.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef AFC_H
#define AFC_H

#include <Material/Material1D/Material1D.h>

struct DataAFC {
    const double elastic_modulus;
    const double t_yield_stress;
    const double t_hardening;
    const double t_unloading;
    const double c_yield_stress;
    const double c_hardening;
    const double c_unloading;
    const double degrade;

    const double t_yield_strain = t_yield_stress / elastic_modulus;
    const double c_yield_strain = c_yield_stress / elastic_modulus;
};

class AFC final : protected DataAFC, public Material1D {
    static podarray<double> compute_transition(double, double, double, double, double, double, double);
    void compute_degradation(double, double);

public:
    AFC(
        unsigned, // unique tag
        double,   // elastic modulus
        double,   // tension yield stress
        double,   // tension hardening modulus
        double,   // tension unloading modulus
        double,   // compression yield stress
        double,   // compression hardening modulus
        double,   // compression unloading modulus
        double,   // degrade parameter
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
