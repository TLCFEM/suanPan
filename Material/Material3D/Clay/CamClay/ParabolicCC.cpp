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

#include "ParabolicCC.h"

double ParabolicCC::compute_a(const double hardening) const { return a_slope < 0. && fabs(hardening) > limit ? 0. : a_slope * hardening * hardening + a; }

double ParabolicCC::compute_da(const double hardening) const { return a_slope < 0. && fabs(hardening) > limit ? 0. : 2. * a_slope * hardening; }

ParabolicCC::ParabolicCC(const unsigned T, const double E, const double V, const double B, const double M, const double P, const double A, const double K, const double R)
    : DataParabolicCC{A, K}
    , NonlinearCamClay(T, E, V, B, M, P, R) {}

unique_ptr<Material> ParabolicCC::get_copy() { return make_unique<ParabolicCC>(*this); }

void ParabolicCC::print() { sp_info("A 3D Cam-Clay model using parabolic hardening.\n"); }
