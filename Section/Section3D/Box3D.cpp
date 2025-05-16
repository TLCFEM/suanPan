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

#include "Box3D.h"

#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Toolbox/IntegrationPlan.h>

Box3D::Box3D(const unsigned T, const double B, const double H, const double TH, const unsigned M, const unsigned S, const double E1, const double E2)
    : Section3D(T, M, 2. * TH * (B + H), vec{E1, E2})
    , width(B)
    , height(H)
    , thickness(TH)
    , int_pt_num(S > 20 ? 20 : S) {}

Box3D::Box3D(const unsigned T, vec&& D, const unsigned M, const unsigned S, vec&& EC)
    : Section3D(T, M, 2. * D(2) * (D(0) + D(1)), std::move(EC))
    , width(D(0))
    , height(D(1))
    , thickness(D(2))
    , int_pt_num(S > 20 ? 20 : S) {}

int Box3D::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get_material(material_tag);

    access::rw(linear_density) = area * material_proto->get_density();

    const IntegrationPlan plan(1, int_pt_num, IntegrationType::GAUSS);

    const auto net_web = height - thickness;
    const auto net_flange = width + thickness;
    const auto web_area = net_web * thickness;
    const auto flange_area = net_flange * thickness;
    const auto web_middle = .5 * width;
    const auto flange_middle = .5 * height;

    int_pt.clear();
    int_pt.reserve(4 * static_cast<size_t>(int_pt_num));
    for(unsigned I = 0; I < int_pt_num; ++I) {
        int_pt.emplace_back(.5 * plan(I, 0) * net_web, web_middle, .5 * plan(I, 1) * web_area, material_proto->get_copy());
        int_pt.emplace_back(.5 * plan(I, 0) * net_web, -web_middle, .5 * plan(I, 1) * web_area, material_proto->get_copy());
        int_pt.emplace_back(flange_middle, .5 * plan(I, 0) * net_flange, .5 * plan(I, 1) * flange_area, material_proto->get_copy());
        int_pt.emplace_back(-flange_middle, .5 * plan(I, 0) * net_flange, .5 * plan(I, 1) * flange_area, material_proto->get_copy());
    }

    initialize_stiffness();

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Box3D::get_copy() { return std::make_unique<Box3D>(*this); }

void Box3D::print() {
    suanpan_info("A 3D box section.\n");
}
