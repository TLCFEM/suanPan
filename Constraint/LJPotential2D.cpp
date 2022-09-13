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

#include "LJPotential2D.h"

double LJPotential2D::compute_f(const double distance) const {
    if(distance >= space) return 0.;

    const auto pow_term = std::pow(.1 * space / distance, 6.);

    return 24. * alpha * pow_term / distance * (2. * pow_term - 1.);
}

double LJPotential2D::compute_df(const double distance) const {
    if(distance >= space) return 0.;

    const auto pow_term = std::pow(.1 * space / distance, 6.);

    return alpha * pow_term * pow(distance, -2.) * (168. - 624. * pow_term);
}
