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
 * @fn NURBSSurface
 * @brief The NURBSSurface class.
 *
 * @author tlc
 * @date 03/11/2020
 * @version 0.1.0
 * @file NURBSSurface.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef NURBSSURFACE_H
#define NURBSSURFACE_H

#include "BSplineSurface.h"
#include "NURBS.h"

class NURBSSurface : public BSplineSurface, public NURBSBase {
public:
    using BSplineSurface::BSplineSurface;

    [[nodiscard]] vec evaluate_point(double, double, const field<vec>&) const override;
    [[nodiscard]] field<vec> evaluate_point_derivative(double, double, const field<vec>&, sword = -1, sword = -1) const override;

    [[nodiscard]] mat evaluate_shape_function(double, double, const field<vec>&) const override;
    [[nodiscard]] field<mat> evaluate_shape_function_derivative(double, double, const field<vec>&, sword = -1, sword = -1) const override;
};

class NURBSSurface2D final : public NURBSSurface {
public:
    explicit NURBSSurface2D(vec, vec, field<vec>&& = {});
};

class NURBSSurface3D final : public NURBSSurface {
public:
    explicit NURBSSurface3D(vec, vec, field<vec>&& = {});
};

class NURBSSurface4D final : public NURBSSurface {
public:
    explicit NURBSSurface4D(vec, vec, field<vec>&& = {});
};

#endif

//! @}
