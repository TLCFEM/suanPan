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

#include "Cell2D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

// eccentricity should really be location and be stored in integration point
// here we flip the sign and use eccentricity to store it
Cell2D::Cell2D(const unsigned T, const double AR, const unsigned MT, const double EC)
    : Section2D(T, MT, AR, -EC) {}

int Cell2D::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get_material(material_tag);

    access::rw(linear_density) = area * material_proto->get_parameter(ParameterType::DENSITY);

    int_pt.clear();
    int_pt.emplace_back(0., area, material_proto->get_copy());

    initialize_stiffness();

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Cell2D::get_copy() { return make_unique<Cell2D>(*this); }

void Cell2D::print() {
    suanpan_info("A 2D section that represents a small cell.\n");
    for(const auto& I : int_pt) I.s_material->print();
}
