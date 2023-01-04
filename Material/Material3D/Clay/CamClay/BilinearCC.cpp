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

#include "BilinearCC.h"

double BilinearCC::compute_a(const double hardening) const { return (a_slope >= 0. ? hardening >= -a / a_slope : hardening < -a / a_slope) ? a + a_slope * hardening : 0.; }

double BilinearCC::compute_da(const double hardening) const { return (a_slope >= 0. ? hardening >= -a / a_slope : hardening < -a / a_slope) ? a_slope : 0.; }

BilinearCC::BilinearCC(const unsigned T, const double E, const double V, const double B, const double M, const double P, const double A, const double K, const double R)
    : DataBilinearCC{A, K}
    , NonlinearCamClay(T, E, V, B, M, P, R) {}

unique_ptr<Material> BilinearCC::get_copy() { return make_unique<BilinearCC>(*this); }

void BilinearCC::print() {
    suanpan_info("A 3D Cam-Clay model using linear hardening.\n");
}
