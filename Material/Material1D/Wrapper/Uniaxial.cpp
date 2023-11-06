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

#include "Uniaxial.h"

Uniaxial::Uniaxial(const unsigned T, const unsigned BT, const unsigned MI)
    : StressWrapper(T, BT, MI, uvec{0}, uvec{1, 2, 3, 4, 5}, MaterialType::D1) {}

unique_ptr<Material> Uniaxial::get_copy() { return make_unique<Uniaxial>(*this); }

void Uniaxial::print() {
    suanpan_info("A uniaxial wrapper.\n");
    suanpan_info("Strain:", current_strain);
    suanpan_info("Stress:", current_stress);
    if(base) base->print();
}
