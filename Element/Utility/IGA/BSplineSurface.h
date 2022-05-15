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
 * @fn BSplineSurface
 * @brief The BSplineSurface class.
 *
 * @author tlc
 * @date 03/11/2020
 * @version 0.1.0
 * @file BSplineSurface.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef BSPLINESURFACE_H
#define BSPLINESURFACE_H

#include "BSpline.h"

class BSplineSurface {
protected:
    const uword dimension;

    field<vec> net;

    BSpline line_u, line_v;

public:
    explicit BSplineSurface(vec, vec, uword, field<vec>&& = {});
    BSplineSurface(const BSplineSurface&) = default;
    BSplineSurface(BSplineSurface&&) = default;
    BSplineSurface& operator=(const BSplineSurface&) = delete;
    BSplineSurface& operator=(BSplineSurface&&) = delete;
    virtual ~BSplineSurface() = default;

    void set_control_polygon(field<vec>&&);
    void set_control_polygon(const field<vec>&);

    [[nodiscard]] field<uvec> get_all_element_span() const;
    [[nodiscard]] uvec get_number_of_control_points() const;

    [[nodiscard]] vec evaluate_point(double, double) const;
    [[nodiscard]] field<vec> evaluate_point_derivative(double, double, sword = -1, sword = -1) const;

    [[nodiscard]] mat evaluate_shape_function(double, double) const;
    [[nodiscard]] field<mat> evaluate_shape_function_derivative(double, double, sword = -1, sword = -1) const;

    [[nodiscard]] virtual vec evaluate_point(double, double, const field<vec>&) const;
    [[nodiscard]] virtual field<vec> evaluate_point_derivative(double, double, const field<vec>&, sword = -1, sword = -1) const;

    [[nodiscard]] virtual mat evaluate_shape_function(double, double, const field<vec>&) const;
    [[nodiscard]] virtual field<mat> evaluate_shape_function_derivative(double, double, const field<vec>&, sword = -1, sword = -1) const;
};

class BSplineSurface2D final : public BSplineSurface {
public:
    explicit BSplineSurface2D(vec, vec, field<vec>&& = {});
};

class BSplineSurface3D final : public BSplineSurface {
public:
    explicit BSplineSurface3D(vec, vec, field<vec>&& = {});
};

class BSplineSurface4D final : public BSplineSurface {
public:
    explicit BSplineSurface4D(vec, vec, field<vec>&& = {});
};

#endif

//! @}
