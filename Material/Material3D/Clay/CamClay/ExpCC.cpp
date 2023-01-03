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

#include "ExpCC.h"

double ExpCC::compute_a(const double hardening) const { return a0 * exp((1. + e0) * (1. - hardening) / (lambda - kappa * hardening)); }

double ExpCC::compute_da(const double hardening) const { return compute_a(hardening) * factor * pow(lambda - hardening * kappa, -2.); }

ExpCC::ExpCC(const unsigned T, const double E, const double V, const double B, const double M, const double P, const double A, const double VR, const double LAMBDA, const double KAPPA, const double R)
    : DataExpCC{A, VR, LAMBDA, KAPPA}
    , NonlinearCamClay(T, E, V, B, M, P, R) {}

unique_ptr<Material> ExpCC::get_copy() { return make_unique<ExpCC>(*this); }

void ExpCC::print() { sp_info("A 3D Cam-Clay model using exponential hardening.\n"); }
