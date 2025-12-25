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
 * @class Balloon
 * @brief A Balloon material class.
 * @author tlc
 * @date 23/12/2025
 * @version 0.1.0
 * @file Balloon.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef BALLOON_H
#define BALLOON_H

#include "BalloonUtil.h"

#include <Material/Material3D/Material3D.h>

struct DataBalloon {
    const double elastic;   // elastic modulus
    const double poisson;   // poisson's ratio
    const double kr;        // plastic strain split ratio
    const unsigned zr_size; // memory size
    const BalloonBuffer::Type zr_type;

    const BalloonBound bound_u, bound_fm, bound_fc, bound_am, bound_ac;

    const std::vector<BalloonSaturation> bfc, bac, bna, bnd;
};

class Balloon final : protected DataBalloon, protected BalloonBase, public Material3D {
    static constexpr unsigned max_iteration = 20u;
    inline static const double root_two_third = std::sqrt(2. / 3.);
    static const mat unit_dev_tensor;

    BalloonBuffer current_zr{zr_size}, trial_zr{zr_size};

    const double double_shear = elastic / (1. + poisson); // double shear modulus

    [[nodiscard]] auto compute_isotropic_bound(double, double, double);
    [[nodiscard]] auto compute_kinematic_bound(double, double, double);

public:
    Balloon(
        unsigned,      // tag
        DataBalloon&&, // data
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
