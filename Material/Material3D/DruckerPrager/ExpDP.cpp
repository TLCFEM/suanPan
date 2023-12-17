/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "ExpDP.h"

double ExpDP::compute_c(const double p_strain) const {
    const auto b_strain = -b * p_strain;
    return cohesion * ((1. + a) * exp(b_strain) - a * exp(2. * b_strain));
}

double ExpDP::compute_dc(const double p_strain) const {
    const auto b_strain = -b * p_strain;
    return cohesion * b * (2. * a * exp(2. * b_strain) - (1. + a) * exp(b_strain));
}

ExpDP::ExpDP(const unsigned T, const double E, const double V, const double ETAY, const double ETAF, const double XI, const double CO, const double PA, const double PB, const double R)
    : DataExpDP({PA <= 1. ? CO : CO * 4. * PA / pow(1. + PA, 2.), PA, PB})
    , NonlinearDruckerPrager(T, E, V, ETAY, ETAF, XI, R) {}

unique_ptr<Material> ExpDP::get_copy() { return make_unique<ExpDP>(*this); }

void ExpDP::print() {
    suanpan_info("A 3D exponential hardening model using Drucker-Prager yielding criterion.\n");
}
