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

#include "Bar3D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

Bar3D::Bar3D(const unsigned T, const double AR, const unsigned MT, const double EA, const double EB)
    : Section3D(T, MT, AR, vec{EA, EB}) {}

int Bar3D::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get_material(material_tag);

    access::rw(linear_density) = area * material_proto->get_parameter(ParameterType::DENSITY);

    int_pt.clear();
    int_pt.emplace_back(0., 0., area, material_proto->get_copy());

    initial_stiffness.set_size(3, 3);

    const auto ea = int_pt.back().s_material->get_initial_stiffness().at(0) * area;
    const auto &arm_y = eccentricity(0), &arm_z = eccentricity(1);
    initial_stiffness(0, 0) = ea;
    initial_stiffness(1, 1) = ea * arm_y * arm_y;
    initial_stiffness(2, 2) = ea * arm_z * arm_z;
    initial_stiffness(1, 0) = initial_stiffness(0, 1) = ea * arm_y;
    initial_stiffness(2, 0) = initial_stiffness(0, 2) = -ea * arm_z;
    initial_stiffness(2, 1) = initial_stiffness(1, 2) = -ea * arm_y * arm_z;

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Bar3D::get_copy() { return make_unique<Bar3D>(*this); }

void Bar3D::print() {
    suanpan_info("A 3D section that represents for example rebar in RC section.\n");
    for(const auto& I : int_pt) I.s_material->print();
}
