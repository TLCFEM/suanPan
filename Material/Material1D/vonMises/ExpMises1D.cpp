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

#include "ExpMises1D.h"

double ExpMises1D::compute_k(const double p_strain) const { return yield_stress + a - a * exp(-b * p_strain) + c * elastic_modulus * p_strain; }

double ExpMises1D::compute_dk(const double p_strain) const { return a * b * exp(-b * p_strain) + c * elastic_modulus; }

double ExpMises1D::compute_h(const double) const { return 0.; }

double ExpMises1D::compute_dh(const double) const { return 0.; }

ExpMises1D::ExpMises1D(const unsigned T, const double E, const double Y, const double A, const double B, const double C, const double R)
    : DataExpMises1D{fabs(Y), fabs(A), fabs(B), fabs(C)}
    , NonlinearMises1D(T, E, R) {}

unique_ptr<Material> ExpMises1D::get_copy() { return std::make_unique<ExpMises1D>(*this); }
