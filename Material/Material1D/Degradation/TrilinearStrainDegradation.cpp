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

#include "TrilinearStrainDegradation.h"

vec TrilinearStrainDegradation::compute_positive_degradation(const double t_strain) const {
    vec damage(2);

    if(const auto abs_e = fabs(t_strain); abs_e > e_strain) {
        damage(0) = e_damage;
        damage(1) = 0.;
    }
    else if(abs_e > s_strain) {
        damage(0) = 1. + slope * (abs_e - s_strain);
        damage(1) = slope;
    }
    else {
        damage(0) = 1.;
        damage(1) = 0.;
    }

    return damage;
}

vec TrilinearStrainDegradation::compute_negative_degradation(const double t_strain) const {
    auto response = compute_positive_degradation(fabs(t_strain));
    response(1) = -response(1);
    return response;
}

TrilinearStrainDegradation::TrilinearStrainDegradation(const unsigned T, const unsigned MT, const double SE, const double EE, const double ED)
    : DataTrilinearStrainDegradation{fabs(SE), fabs(EE), fabs(ED)}
    , StrainDegradation(T, MT) {}

unique_ptr<Material> TrilinearStrainDegradation::get_copy() { return std::make_unique<TrilinearStrainDegradation>(*this); }
