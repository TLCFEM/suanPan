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

#include "BilinearMises1D.h"

double BilinearMises1D::compute_k(const double p_strain) const { return yield_stress + p_strain * isotropic_modulus; }

double BilinearMises1D::compute_dk(const double) const { return isotropic_modulus; }

double BilinearMises1D::compute_h(const double) const { return 0.; }

double BilinearMises1D::compute_dh(const double) const { return 0.; }

BilinearMises1D::BilinearMises1D(const unsigned T, const double E, const double Y, const double H, const double R)
    : DataBilinearMises1D{fabs(Y), fabs(E) * H / (1. - H)}
    , NonlinearMises1D(T, E, R) {}

unique_ptr<Material> BilinearMises1D::get_copy() { return make_unique<BilinearMises1D>(*this); }
