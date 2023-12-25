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

#include "Rectangle3D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Toolbox/IntegrationPlan.h>

Rectangle3D::Rectangle3D(const unsigned T, const double B, const double H, const unsigned M, const unsigned S, const double E1, const double E2)
    : Section3D(T, M, B * H, vec{E1, E2})
    , width(B)
    , height(H)
    , int_pt_num(S > 20 ? 20 : S) {}

int Rectangle3D::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get_material(material_tag);

    access::rw(linear_density) = area * material_proto->get_density();

    const IntegrationPlan plan_y(1, int_pt_num, IntegrationType::LOBATTO);
    const IntegrationPlan plan_z(1, int_pt_num, IntegrationType::LOBATTO);

    int_pt.clear();
    int_pt.reserve(static_cast<size_t>(int_pt_num) * static_cast<size_t>(int_pt_num));
    for(unsigned I = 0; I < int_pt_num; ++I) for(unsigned J = 0; J < int_pt_num; ++J) int_pt.emplace_back(.5 * height * plan_y(I, 0), .5 * width * plan_z(J, 0), .25 * plan_y(I, 1) * plan_z(J, 1) * area, material_proto->get_copy());

    initialize_stiffness();

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Rectangle3D::get_copy() { return make_unique<Rectangle3D>(*this); }

void Rectangle3D::print() {
    suanpan_info("A 3D rectangular section.\n");
}
