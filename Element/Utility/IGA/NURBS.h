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
 * @fn NURBS
 * @brief The NURBS class.
 *
 * @author tlc
 * @date 03/11/2020
 * @version 0.1.0
 * @file NURBS.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef NURBS_H
#define NURBS_H

#include "BSpline.h"

class NURBSBase {
protected:
    mat binomial_mat;

    void initialize_binomial(sword) const;
};

class NURBS : public BSpline, public NURBSBase {
public:
    using BSpline::BSpline;

    [[nodiscard]] vec evaluate_point(double, const field<vec>&) const override;
    [[nodiscard]] field<vec> evaluate_point_derivative(double, const field<vec>&, sword = -1) const override;

    [[nodiscard]] vec evaluate_shape_function(double, const field<vec>&) const override;
    [[nodiscard]] field<vec> evaluate_shape_function_derivative(double, const field<vec>&, sword = -1) const override;
};

class NURBSCurve2D final : public NURBS {
public:
    explicit NURBSCurve2D(vec, field<vec>&& = {});
};

class NURBSCurve3D final : public NURBS {
public:
    explicit NURBSCurve3D(vec, field<vec>&& = {});
};

class NURBSCurve4D final : public NURBS {
public:
    explicit NURBSCurve4D(vec, field<vec>&& = {});
};

#endif

//! @}
