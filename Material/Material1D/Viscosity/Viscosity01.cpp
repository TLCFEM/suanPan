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

#include "Viscosity01.h"

double Viscosity01::compute_du(double, double) const { return 0.; }

double Viscosity01::compute_dv(double, double) const { return 0.; }

double Viscosity01::compute_damping_coefficient(double, double) const { return damping; }

Viscosity01::Viscosity01(const unsigned T, const double A, const double C, const double L)
    : DataViscosity01{fabs(C)}
    , NonlinearViscosity(T, A, L) {}

unique_ptr<Material> Viscosity01::get_copy() { return make_unique<Viscosity01>(*this); }
