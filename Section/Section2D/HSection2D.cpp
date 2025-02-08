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

#include "HSection2D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Toolbox/IntegrationPlan.h>

HSection2D::HSection2D(const unsigned T, const double TFW, const double TFT, const double BFW, const double BFT, const double WH, const double WT, const unsigned MT, const unsigned IP, const double EC)
    : Section2D(T, MT, TFW * TFT + BFW * BFT + WH * WT, EC)
    , left_flange_height(TFW)
    , left_flange_thickness(TFT)
    , right_flange_height(BFW)
    , right_flange_thickness(BFT)
    , web_width(WH)
    , web_thickness(WT)
    , int_pt_num(IP > 20 ? 20 : IP) {}

int HSection2D::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get_material(material_tag);

    access::rw(linear_density) = mat_proto->get_density() * area;

    const auto web_area = web_width * web_thickness;
    const auto l_flange_area = left_flange_height * left_flange_thickness;
    const auto r_flange_area = right_flange_height * right_flange_thickness;

    const IntegrationPlan plan_flange(1, int_pt_num, IntegrationType::GAUSS);
    const IntegrationPlan plan_web(1, 2, IntegrationType::GAUSS);

    int_pt.clear();
    int_pt.reserve(2llu * int_pt_num + 2);
    int_pt.emplace_back(.5 * plan_web(0, 0) * web_thickness, .5 * plan_web(0, 1) * web_area, mat_proto->get_copy());
    int_pt.emplace_back(.5 * plan_web(1, 0) * web_thickness, .5 * plan_web(1, 1) * web_area, mat_proto->get_copy());
    for(unsigned I = 0; I < int_pt_num; ++I) {
        int_pt.emplace_back(.5 * plan_flange(I, 0) * left_flange_height, .5 * plan_flange(I, 1) * l_flange_area, mat_proto->get_copy());
        int_pt.emplace_back(.5 * plan_flange(I, 0) * right_flange_height, .5 * plan_flange(I, 1) * r_flange_area, mat_proto->get_copy());
    }

    initialize_stiffness();

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> HSection2D::get_copy() { return make_unique<HSection2D>(*this); }

void HSection2D::print() {
    suanpan_info("A 2D H-shape section.\n");
    auto J = 1;
    for(const auto& I : int_pt) {
        suanpan_info("IP {}: {:.4E}.\n", J++, I.coor);
        I.s_material->print();
    }
}
