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
 * @fn BSpline
 * @brief The BSpline class.
 *
 * @author tlc
 * @date 03/11/2020
 * @version 0.1.0
 * @file BSpline.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef BSPLINE_H
#define BSPLINE_H

#include <suanPan.h>

struct IGA {
    static void convert_to_weighted(mat&);
    static void convert_to_weighted(field<vec>&);
    [[nodiscard]] static uword compute_order(const vec&);
    [[nodiscard]] static uword compute_number_of_elements(const vec&);
    [[nodiscard]] static uword compute_number_of_control_points(const vec&);
    [[nodiscard]] static uvec compute_all_element_span(const vec&);
};

class BSpline {
protected:
    const uword dimension;

    const vec knot;

    field<vec> net;

    const uword order = IGA::compute_order(knot);

public:
    BSpline(vec, uword, field<vec>&& = {});
    BSpline(const BSpline&) = default;
    BSpline(BSpline&&) = default;
    BSpline& operator=(const BSpline&) = delete;
    BSpline& operator=(BSpline&&) = delete;
    virtual ~BSpline() = default;

    void set_control_polygon(field<vec>&&);
    void set_control_polygon(const field<vec>&);

    [[nodiscard]] const vec& get_knot() const;
    [[nodiscard]] uword get_order() const;
    [[nodiscard]] uword get_number_of_control_points() const;
    [[nodiscard]] uvec get_all_element_span() const;

    [[nodiscard]] uword evaluate_span(double) const;

    [[nodiscard]] vec evaluate_basis(double, sword = -1) const;
    [[nodiscard]] mat evaluate_basis_derivative(double, sword = -1, sword = -1) const;

    [[nodiscard]] vec evaluate_point(double) const;
    [[nodiscard]] field<vec> evaluate_point_derivative(double, sword = -1) const;

    [[nodiscard]] vec evaluate_shape_function(double) const;
    [[nodiscard]] field<vec> evaluate_shape_function_derivative(double, sword = -1) const;

    [[nodiscard]] virtual vec evaluate_point(double, const field<vec>&) const;
    [[nodiscard]] virtual field<vec> evaluate_point_derivative(double, const field<vec>&, sword = -1) const;

    [[nodiscard]] virtual vec evaluate_shape_function(double, const field<vec>&) const;
    [[nodiscard]] virtual field<vec> evaluate_shape_function_derivative(double, const field<vec>&, sword = -1) const;
};

class BSplineCurve2D final : public BSpline {
public:
    explicit BSplineCurve2D(vec, field<vec>&& = {});
};

class BSplineCurve3D final : public BSpline {
public:
    explicit BSplineCurve3D(vec, field<vec>&& = {});
};

class BSplineCurve4D final : public BSpline {
public:
    explicit BSplineCurve4D(vec, field<vec>&& = {});
};

#endif

//! @}
