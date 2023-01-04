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

#include "BilinearDP.h"

double BilinearDP::compute_c(const double p_strain) const { return cohesion_slope >= 0. || p_strain <= critical ? cohesion + cohesion_slope * p_strain : 0.; }

double BilinearDP::compute_dc(const double p_strain) const { return cohesion_slope >= 0. || p_strain <= critical ? cohesion_slope : 0.; }

BilinearDP::BilinearDP(const unsigned T, const double E, const double V, const double ETAY, const double ETAF, const double XI, const double CO, const double CS, const double R)
    : DataBilinearDP{CO, CS}
    , NonlinearDruckerPrager(T, E, V, ETAY, ETAF, XI, R) {}

unique_ptr<Material> BilinearDP::get_copy() { return make_unique<BilinearDP>(*this); }

void BilinearDP::print() {
    suanpan_info("A 3D nonlinear model using Drucker-Prager yielding criterion with linear cohesion.\n");
}
