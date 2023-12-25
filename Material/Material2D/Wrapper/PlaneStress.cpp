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

#include "PlaneStress.h"
#include <Domain/DomainBase.h>

PlaneStress::PlaneStress(const unsigned T, const unsigned BT, const unsigned MI)
    : StressWrapper(T, BT, MI, uvec{0, 1, 3}, uvec{2, 4, 5}, MaterialType::D2) { access::rw(plane_type) = PlaneType::S; }

unique_ptr<Material> PlaneStress::get_copy() { return make_unique<PlaneStress>(*this); }

void PlaneStress::print() {
    suanpan_info("A plane stress wrapper.\n");
    suanpan_info("Strain:", current_strain);
    suanpan_info("Stress:", current_stress);
    if(base) base->print();
}
