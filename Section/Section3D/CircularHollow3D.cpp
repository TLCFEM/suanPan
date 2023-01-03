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

#include "CircularHollow3D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Toolbox/IntegrationPlan.h>

CircularHollow3D::CircularHollow3D(const unsigned T, const double R, const double TH, const unsigned M, const unsigned S, vec&& EC)
    : Section3D(T, M, (R * R - (R - TH) * (R - TH)) * datum::pi, std::forward<vec>(EC))
    , radius(R)
    , thickness(TH)
    , int_pt_num(S > 20 ? 20 : S) {}

CircularHollow3D::CircularHollow3D(const unsigned T, vec&& D, const unsigned M, const unsigned S, vec&& EC)
    : Section3D(T, M, (D(0) * D(0) - (D(0) - D(1)) * (D(0) - D(1))) * datum::pi, std::forward<vec>(EC))
    , radius(D(0))
    , thickness(D(1))
    , int_pt_num(S > 20 ? 20 : S) {}

int CircularHollow3D::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get_material(material_tag);

    access::rw(linear_density) = area * material_proto->get_parameter(ParameterType::DENSITY);

    const IntegrationPlan plan(1, int_pt_num, IntegrationType::GAUSS);

    const auto m_radius = radius - .5 * thickness;

    int_pt.clear();
    int_pt.reserve(2llu * int_pt_num);
    initial_stiffness.zeros(3, 3);
    for(unsigned I = 0; I < int_pt_num; ++I) {
        const auto t_angle = .5 * plan(I, 0) * datum::pi;
        int_pt.emplace_back(cos(t_angle) * m_radius, sin(t_angle) * m_radius, .25 * plan(I, 1) * area, material_proto->get_copy());
        int_pt.emplace_back(-cos(t_angle) * m_radius, -sin(t_angle) * m_radius, .25 * plan(I, 1) * area, material_proto->get_copy());
    }
    for(const auto& I : int_pt) {
        const auto arm_y = eccentricity(0) - I.coor_y;
        const auto arm_z = I.coor_z - eccentricity(1);
        const auto tmp_a = I.s_material->get_initial_stiffness().at(0) * I.weight;
        initial_stiffness(0, 0) += tmp_a;
        initial_stiffness(0, 1) += tmp_a * arm_y;
        initial_stiffness(0, 2) += tmp_a * arm_z;
        initial_stiffness(1, 1) += tmp_a * arm_y * arm_y;
        initial_stiffness(1, 2) += tmp_a * arm_y * arm_z;
        initial_stiffness(2, 2) += tmp_a * arm_z * arm_z;
    }
    initial_stiffness(1, 0) = initial_stiffness(0, 1);
    initial_stiffness(2, 0) = initial_stiffness(0, 2);
    initial_stiffness(2, 1) = initial_stiffness(1, 2);

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> CircularHollow3D::get_copy() { return make_unique<CircularHollow3D>(*this); }

void CircularHollow3D::print() { sp_info("A 3D circular hollow section.\n"); }
