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

#include "Dhakal.h"

podarray<double> Dhakal::compute_degradation(const double t_strain) const {
    podarray<double> damage(2);

    if(t_strain < final_strain) {
        damage(0) = .2;
        damage(1) = 0.;
    }
    else if(t_strain < inter_strain) {
        damage(1) = -.02 / yield_strain;
        damage(0) = (t_strain - inter_strain) * damage(1) + inter_factor;
    }
    else if(t_strain < yield_strain) {
        damage(1) = slope;
        damage(0) = (t_strain - yield_strain) * damage(1) + 1.;
    }
    else {
        damage(0) = 1.;
        damage(1) = 0.;
    }

    return damage;
}

Dhakal::Dhakal(const unsigned T, const unsigned MT, const double EY, const double PP)
    : DataDhakal{-fabs(EY), -fabs(EY) * std::max(55. - 2.3 * fabs(PP), 7.), std::min(1., std::max(1.1 - .016 * fabs(PP), .2))}
    , Degradation(T, MT) {}

unique_ptr<Material> Dhakal::get_copy() { return make_unique<Dhakal>(*this); }
