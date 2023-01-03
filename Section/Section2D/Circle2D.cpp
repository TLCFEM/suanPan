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

#include "Circle2D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Toolbox/IntegrationPlan.h>

Circle2D::Circle2D(const unsigned T, const double R, const unsigned M, const unsigned S, const double EC)
    : Section2D(T, M, R * R * datum::pi, EC)
    , radius(R)
    , int_pt_num(S > 20 ? 20 : S) {}

int Circle2D::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get_material(material_tag);

    access::rw(linear_density) = area * material_proto->get_parameter(ParameterType::DENSITY);

    const IntegrationPlan plan(1, int_pt_num, IntegrationType::GAUSS);

    int_pt.clear();
    int_pt.reserve(int_pt_num);
    initial_stiffness.zeros(2, 2);
    for(unsigned I = 0; I < int_pt_num; ++I) {
        int_pt.emplace_back(radius * plan(I, 0), 2. * radius * radius * sqrt(1. - plan(I, 0) * plan(I, 0)) * plan(I, 1), material_proto->get_copy());
        auto tmp_a = int_pt[I].s_material->get_initial_stiffness().at(0) * int_pt[I].weight;
        const auto arm = eccentricity(0) - int_pt[I].coor;
        initial_stiffness(0, 0) += tmp_a;
        initial_stiffness(0, 1) += tmp_a *= arm;
        initial_stiffness(1, 1) += tmp_a *= arm;
    }
    initial_stiffness(1, 0) = initial_stiffness(0, 1);

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Circle2D::get_copy() { return make_unique<Circle2D>(*this); }

void Circle2D::print() { sp_info("A 2D circular section.\n"); }
