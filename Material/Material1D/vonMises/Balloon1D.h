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
 * @class Balloon1D
 * @brief A Balloon1D material class.
 * @author tlc
 * @date 04/11/2025
 * @version 0.1.0
 * @file Balloon1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef BALLOON1D_H
#define BALLOON1D_H

#include <Material/Material1D/Material1D.h>
#include <Material/Material3D/vonMises/BalloonUtil.h>

struct DataBalloon1D {
    const double elastic;   // elastic modulus
    const double kr;        // plastic strain split ratio
    const unsigned zr_size; // memory size
    const BalloonBuffer::Type zr_type;

    const BalloonBound bound_u, bound_fm, bound_fc, bound_am, bound_ac;

    const std::vector<BalloonSaturation> bfc, bac, bna, bnd;
};

class Balloon1D final : protected DataBalloon1D, protected BalloonBase, public Material1D {
    static constexpr unsigned max_iteration = 20u;

    BalloonBuffer current_zr{zr_size}, trial_zr{zr_size};

    [[nodiscard]] double initial_check(double);

    [[nodiscard]] auto compute_isotropic_bound(double, double, double);
    [[nodiscard]] auto compute_kinematic_bound(double, double, double);

public:
    Balloon1D(
        unsigned,        // tag
        DataBalloon1D&&, // data
        double = 0.      // density
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
