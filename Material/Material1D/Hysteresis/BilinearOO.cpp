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

#include "BilinearOO.h"

podarray<double> BilinearOO::compute_tension_backbone(const double strain) const {
    podarray<double> response(2);

    response(0) = strain > t_strain ? t_stress + (response(1) = t_hardening) * (strain - t_strain) : (response(1) = elastic_modulus) * strain;

    if(response(0) < 0.) {
        response(0) = 0.;
        response(1) = 1E-10;
    }

    return response;
}

podarray<double> BilinearOO::compute_compression_backbone(const double strain) const {
    podarray<double> response(2);

    response(0) = strain < c_strain ? c_stress + (response(1) = c_hardening) * (strain - c_strain) : (response(1) = elastic_modulus) * strain;

    if(response(0) > 0.) {
        response(0) = 0.;
        response(1) = 1E-10;
    }

    return response;
}

BilinearOO::BilinearOO(const int T, const double E, const double TEA, const double TH, const double CEA, const double CH, const double R)
    : DataBilinearOO{fabs(E), fabs(TEA), TH * fabs(E), -fabs(CEA), CH * fabs(E)}
    , OriginOriented(T, R) {}

int BilinearOO::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> BilinearOO::get_copy() { return make_unique<BilinearOO>(*this); }

void BilinearOO::print() {
    suanpan_info("A bilinear origin oriented hysteresis model.\n");
    OriginOriented::print();
}
