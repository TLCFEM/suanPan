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

#include "Cell3DOS.h"

#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

Cell3DOS::Cell3DOS(const unsigned T, const double AR, const double OM, const double PY, const double PZ, const unsigned MT, const double EA, const double EB)
    : SectionOS3D(T, MT, AR, vec{EA, EB})
    , omega(OM)
    , py(PY)
    , pz(PZ) {}

int Cell3DOS::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get_material(material_tag);

    access::rw(linear_density) = area * material_proto->get_density();

    int_pt.clear();
    int_pt.emplace_back(omega, py, pz, area, material_proto->get_copy());

    const auto& arm_y = eccentricity(0);
    const auto& arm_z = eccentricity(1);

    const auto os_size = static_cast<unsigned>(section_type);

    sp_mat de(3, os_size);
    de(0, 0) = 1.;
    de(0, 3) = -arm_y;
    de(0, 4) = -arm_z;
    de(0, 7) = omega;
    de(1, 6) = py - arm_z;
    de(2, 6) = pz + arm_y;

    trial_stiffness = current_stiffness = initial_stiffness = area * de.t() * int_pt.back().s_material->get_initial_stiffness() * de;

    trial_geometry = current_geometry = initial_geometry.zeros(os_size, os_size);

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Cell3DOS::get_copy() { return std::make_unique<Cell3DOS>(*this); }

void Cell3DOS::print() {
    suanpan_info("A 3D open section cell.\n");
    for(const auto& I : int_pt) I.s_material->print();
}
