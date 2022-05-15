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
 * @fn BSplineVolume
 * @brief The BSplineVolume class.
 *
 * @author tlc
 * @date 03/11/2020
 * @version 0.1.0
 * @file BSplineVolume.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef BSPLINEVOLUME_H
#define BSPLINEVOLUME_H

#include "BSpline.h"

class BSplineVolume {
protected:
    const uword dimension;

    field<vec> net;

    BSpline line_u, line_v, line_w;

public:
    explicit BSplineVolume(vec, vec, vec, uword, field<vec>&& = {});
    BSplineVolume(const BSplineVolume&) = default;
    BSplineVolume(BSplineVolume&&) = default;
    BSplineVolume& operator=(const BSplineVolume&) = delete;
    BSplineVolume& operator=(BSplineVolume&&) = delete;
    virtual ~BSplineVolume() = default;

    void set_control_polygon(field<vec>&&);
    void set_control_polygon(const field<vec>&);

    [[nodiscard]] field<uvec> get_all_element_span() const;
    [[nodiscard]] uvec get_number_of_control_points() const;

    [[nodiscard]] vec evaluate_point(double, double, double) const;
    [[nodiscard]] field<vec> evaluate_point_derivative(double, double, double, sword = -1) const;

    [[nodiscard]] cube evaluate_shape_function(double, double, double) const;
    [[nodiscard]] field<cube> evaluate_shape_function_derivative(double, double, double, sword = -1, sword = -1, sword = -1) const;

    [[nodiscard]] virtual vec evaluate_point(double, double, double, const field<vec>&) const;
    [[nodiscard]] virtual field<vec> evaluate_point_derivative(double, double, double, const field<vec>&, sword = -1) const;

    [[nodiscard]] virtual cube evaluate_shape_function(double, double, double, const field<vec>&) const;
    [[nodiscard]] virtual field<cube> evaluate_shape_function_derivative(double, double, double, const field<vec>&, sword = -1, sword = -1, sword = -1) const;
};

class BSplineVolume3D final : public BSplineVolume {
public:
    explicit BSplineVolume3D(vec, vec, vec, field<vec>&& = {});
};

class BSplineVolume4D final : public BSplineVolume {
public:
    explicit BSplineVolume4D(vec, vec, vec, field<vec>&& = {});
};

#endif

//! @}
