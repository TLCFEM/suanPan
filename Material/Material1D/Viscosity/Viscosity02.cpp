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

#include "Viscosity02.h"

double Viscosity02::compute_du(const double strain, const double strain_rate) const { return (damping_b + damping_d * atan(gap_b * strain_rate)) * gap_a / (1. + pow(gap_a * strain, 2.)); }

double Viscosity02::compute_dv(const double strain, const double strain_rate) const { return (damping_c + damping_d * atan(gap_a * strain)) * gap_b / (1. + pow(gap_b * strain_rate, 2.)); }

double Viscosity02::compute_damping_coefficient(const double strain, const double strain_rate) const {
    const auto factor_a = atan(gap_a * strain);
    const auto factor_b = atan(gap_b * strain_rate);
    return damping_a + damping_b * factor_a + damping_c * factor_b + damping_d * factor_a * factor_b;
}

Viscosity02::Viscosity02(const unsigned T, const double A, const double CA, const double CB, const double CC, const double CD, const double GA, const double GB, const double L)
    : DataViscosity02{.25 * (fabs(CA) + fabs(CB) + fabs(CC) + fabs(CD)), .5 / datum::pi * (fabs(CA) - fabs(CB) - fabs(CC) + fabs(CD)), .5 / datum::pi * (fabs(CA) + fabs(CB) - fabs(CC) - fabs(CD)), (fabs(CA) - fabs(CB) + fabs(CC) - fabs(CD)) / datum::pi / datum::pi, fabs(GA), fabs(GB)}
    , NonlinearViscosity(T, fabs(A), fabs(L)) {}

unique_ptr<Material> Viscosity02::get_copy() { return make_unique<Viscosity02>(*this); }
