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

#include "Circle3D.h"

#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Toolbox/IntegrationPlan.h>

Circle3D::Circle3D(const unsigned T, const double R, const unsigned M, const unsigned S, vec&& EC)
    : Section3D(T, M, R * R * datum::pi, std::move(EC))
    , radius(R)
    , int_pt_num(S > 20 ? 20 : S) {}

int Circle3D::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get_material(material_tag);

    access::rw(linear_density) = area * material_proto->get_density();

    const IntegrationPlan plan(2, int_pt_num, IntegrationType::GAUSS);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const auto t_angle = (plan(I, 0) + 1.) * datum::pi;
        const auto t_radius = .5 * radius * (plan(I, 1) + 1.);
        int_pt.emplace_back(cos(t_angle) * t_radius, sin(t_angle) * t_radius, .5 * plan(I, 2) * t_radius * area, material_proto->get_copy());
    }

    initialize_stiffness();

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Circle3D::get_copy() { return std::make_unique<Circle3D>(*this); }

void Circle3D::print() {
    suanpan_info("A 3D circular section.\n");
}
