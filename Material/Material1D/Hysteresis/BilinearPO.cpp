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

#include "BilinearPO.h"

pod2 BilinearPO::compute_compression_initial_reverse() const {
    pod2 response;

    response[1] = (response[0] = c_strain) * initial_stiffness(0);

    return response;
}

pod2 BilinearPO::compute_tension_initial_reverse() const {
    pod2 response;

    response[1] = (response[0] = t_strain) * initial_stiffness(0);

    return response;
}

pod2 BilinearPO::compute_tension_backbone(const double strain) const {
    pod2 response;

    response[0] = strain > t_strain ? t_stress + (response[1] = t_hardening) * (strain - t_strain) : (response[1] = elastic_modulus) * strain;

    if(response[0] < 0.) {
        response[0] = 0.;
        response[1] = 1E-10;
    }

    return response;
}

pod2 BilinearPO::compute_compression_backbone(const double strain) const {
    pod2 response;

    response[0] = strain < c_strain ? c_stress + (response[1] = c_hardening) * (strain - c_strain) : (response[1] = elastic_modulus) * strain;

    if(response[0] > 0.) {
        response[0] = 0.;
        response[1] = 1E-10;
    }

    return response;
}

BilinearPO::BilinearPO(const int T, const double E, const double TEA, const double TH, const double CEA, const double CH, const double R)
    : DataBilinearPO{fabs(E), fabs(TEA), TH * fabs(E), -fabs(CEA), CH * fabs(E)}
    , PeakOriented(T, R) {}

int BilinearPO::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> BilinearPO::get_copy() { return make_unique<BilinearPO>(*this); }

void BilinearPO::print() {
    suanpan_info("A bilinear peak oriented hysteresis model.\n");
    PeakOriented::print();
}
