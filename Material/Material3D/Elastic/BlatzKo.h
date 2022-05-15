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
 * @class BlatzKo
 * @brief The BlatzKo class.
 *
 * algorithm verified at 24 April 2019 by tlc
 *
 * @author tlc
 * @date 24/04/2019
 * @version 1.0.0
 * @file BlatzKo.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef BLATZKO_H
#define BLATZKO_H

#include <Material/Material3D/Material3D.h>

struct DataBlatzKo {
    const double elastic_modulus;
    const double poissons_ratio;
    const double shear_modulus;
    const double half_beta_two;
};

class BlatzKo final : DataBlatzKo, public Material3D {
    static const vec weight;

public:
    BlatzKo(unsigned,   // tag
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
