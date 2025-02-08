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
 * @fn NURBSVolume
 * @brief The NURBSVolume class.
 *
 * @author tlc
 * @date 03/11/2020
 * @version 0.1.0
 * @file NURBSVolume.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef NURBSVOLUME_H
#define NURBSVOLUME_H

#include "BSplineVolume.h"
#include "NURBS.h"

class NURBSVolume : public BSplineVolume, public NURBSBase {
public:
    using BSplineVolume::BSplineVolume;

    [[nodiscard]] vec evaluate_point(double, double, double, const field<vec>&) const override;
    [[nodiscard]] field<vec> evaluate_point_derivative(double, double, double, const field<vec>&, sword = -1) const override;

    [[nodiscard]] cube evaluate_shape_function(double, double, double, const field<vec>&) const override;
    [[nodiscard]] field<cube> evaluate_shape_function_derivative(double, double, double, const field<vec>&, sword = -1, sword = -1, sword = -1) const override;
};

class NURBSVolume3D final : public NURBSVolume {
public:
    explicit NURBSVolume3D(vec, vec, vec, field<vec>&& = {});
};

class NURBSVolume4D final : public NURBSVolume {
public:
    explicit NURBSVolume4D(vec, vec, vec, field<vec>&& = {});
};

#endif

//! @}
