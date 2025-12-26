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
 * @class Subloading1D
 * @brief A Subloading1D material class.
 * @author tlc
 * @date 24/09/2024
 * @version 0.1.0
 * @file Subloading1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef SUBLOADING1D_H
#define SUBLOADING1D_H

#include <Material/Material1D/Material1D.h>
#include <Material/Material3D/vonMises/SubloadingUtil.h>

struct DataSubloading1D {
    const double elastic; // elastic modulus
    const SubloadingBound iso_bound, kin_bound;
    const double u;
    const double cv;
    const double mu;
    const double nv;

    const std::vector<SubloadingSaturation> b, c;
};

class Subloading1D final : protected DataSubloading1D, protected SubloadingBase, public Material1D {
    static constexpr unsigned max_iteration = 20u;

    const double* incre_time = nullptr;

    const bool is_viscous = mu > 0. && nv > 0.;

public:
    Subloading1D(
        unsigned,           // tag
        DataSubloading1D&&, // data
        double = 0.         // density
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
