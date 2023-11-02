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
 * @class Yeoh
 * @brief The Yeoh class.
 *
 * @author tlc
 * @date 25/12/2020
 * @version 1.0.0
 * @file Yeoh.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef YEOH_H
#define YEOH_H

#include <Material/Material3D/Material3D.h>

struct DataYeoh {
    const vec A0;
    const vec A1;
};

class Yeoh final : protected DataYeoh, public Material3D {
    static const vec weight;
    static const vec I1E;

    static constexpr double one_three = 1. / 3.;
    static constexpr double two_three = 2. * one_three;
    static constexpr double four_three = 2. * two_three;
    static constexpr double five_three = 5. * one_three;
    static constexpr double eight_nine = two_three * four_three;

    [[nodiscard]] vec compute_derivative(double, double) const;

public:
    Yeoh(
        unsigned,   // tag
        vec&&,      // constants
        vec&&,      // constants
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
