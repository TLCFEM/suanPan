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
 * @class SurfaceNM2D
 * @brief A SurfaceNM2D class.
 * @author tlc
 * @date 22/06/2022
 * @version 0.1.0
 * @file SurfaceNM2D.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef SURFACENM2D_H
#define SURFACENM2D_H

#include <suanPan.h>

class SurfaceNM2D {
    const mat para_set;
    const double c;

    [[nodiscard]] static double evaluate(double, double, const mat&);
    [[nodiscard]] static vec differentiate(const mat&, uword, uword);

public:
    explicit SurfaceNM2D(
        double, // c
        mat&& = {}
    );

    [[nodiscard]] double compute_sf(const vec&, const vec&) const;
    [[nodiscard]] vec compute_dsf(const vec&, const vec&) const;
    [[nodiscard]] mat compute_ddsf(const vec&, const vec&) const;
};

#endif

//! @}
