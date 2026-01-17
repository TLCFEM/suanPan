/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
 * @class Subloading
 * @brief A Subloading material class.
 * @author tlc
 * @date 26/12/2025
 * @version 0.2.0
 * @file Subloading.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef SUBLOADING_H
#define SUBLOADING_H

#include "SubloadingUtil.h"

#include <Material/Material3D/Material3D.h>

struct DataSubloading {
    const double elastic; // elastic modulus
    const double poisson;
    const SubloadingBound iso_bound, kin_bound;
    const double u;

    const SubloadingSaturation b, c;
};

class Subloading final : protected DataSubloading, protected SubloadingBase, public Material3D {
    static constexpr unsigned max_iteration = 20u;
    inline static const double root_two_third = std::sqrt(2. / 3.);
    static const mat unit_dev_tensor;

    const double double_shear = elastic / (1. + poisson); // double shear modulus

public:
    Subloading(
        unsigned,         // tag
        DataSubloading&&, // data
        double = 0.       // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> unique_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
