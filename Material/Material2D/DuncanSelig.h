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
 * @class DuncanSelig
 * @brief A DuncanSelig material class.
 *
 * For plane strain soil.
 *
 * @author tlc
 * @date 01/04/24
 * @version 1.0.0
 * @file DuncanSelig.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef DUNCANSELIG_H
#define DUNCANSELIG_H

#include <Material/Material2D/Material2D.h>

class DuncanSelig final : public Material2D {
    const double elastic_modulus;
    double shear_modulus = 0.;

public:
    DuncanSelig(
        unsigned,   // tag
        double,     // elastic modulus
        double = 0. // density
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