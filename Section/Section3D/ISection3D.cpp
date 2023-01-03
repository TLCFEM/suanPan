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

#include "ISection3D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Toolbox/IntegrationPlan.h>

ISection3D::ISection3D(const unsigned T, const double TFW, const double TFT, const double BFW, const double BFT, const double WH, const double WT, const unsigned MT, const unsigned IP, vec&& EC)
    : Section3D(T, MT, TFW * TFT + BFW * BFT + WH * WT, std::forward<vec>(EC))
    , top_flange_width(TFW)
    , top_flange_thickness(TFT)
    , bottom_flange_width(BFW)
    , bottom_flange_thickness(BFT)
    , web_height(WH)
    , web_thickness(WT)
    , int_pt_num(IP > 20 ? 20 : IP) {}

ISection3D::ISection3D(const unsigned T, vec&& D, const unsigned MT, const unsigned IP, vec&& EC)
    : Section3D(T, MT, D(0) * D(1) + D(2) * D(3) + D(4) * D(5), std::forward<vec>(EC))
    , top_flange_width(D(0))
    , top_flange_thickness(D(1))
    , bottom_flange_width(D(2))
    , bottom_flange_thickness(D(3))
    , web_height(D(4))
    , web_thickness(D(5))
    , int_pt_num(IP > 20 ? 20 : IP) {}

int ISection3D::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get_material(material_tag);

    access::rw(linear_density) = mat_proto->get_parameter(ParameterType::DENSITY) * area;

    const auto web_area = web_height * web_thickness;
    const auto b_flange_area = bottom_flange_width * bottom_flange_thickness;
    const auto t_flange_area = top_flange_width * top_flange_thickness;

    const IntegrationPlan plan_flange(1, int_pt_num, IntegrationType::GAUSS);
    const auto& plan_web = plan_flange;

    int_pt.clear();
    int_pt.reserve(plan_web.n_rows + 2llu * plan_flange.n_rows);
    for(unsigned I = 0; I < int_pt_num; ++I) int_pt.emplace_back(.5 * plan_web(I, 0) * web_height, 0., .5 * plan_web(I, 1) * web_area, mat_proto->get_copy());
    if(b_flange_area != 0.) for(unsigned I = 0; I < plan_flange.n_rows; ++I) int_pt.emplace_back(.5 * (bottom_flange_thickness + web_height), .5 * plan_flange(I, 0) * bottom_flange_width, .5 * plan_flange(I, 1) * b_flange_area, mat_proto->get_copy());
    if(t_flange_area != 0.) for(unsigned I = 0; I < plan_flange.n_rows; ++I) int_pt.emplace_back(-.5 * (top_flange_thickness + web_height), .5 * plan_flange(I, 0) * top_flange_width, .5 * plan_flange(I, 1) * t_flange_area, mat_proto->get_copy());

    initial_stiffness.zeros(3, 3);
    for(const auto& I : int_pt) {
        const auto tmp_a = I.s_material->get_initial_stiffness().at(0) * I.weight;
        const auto arm_y = eccentricity(0) - I.coor_y;
        const auto arm_z = I.coor_z - eccentricity(1);
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

unique_ptr<Section> ISection3D::get_copy() { return make_unique<ISection3D>(*this); }

void ISection3D::print() {
    sp_info("A 3D I-shape section with following integration points.\n");
    auto J = 1;
    for(const auto& I : int_pt) {
        sp_info("IP {}: {:.4E}, {:.4E}.\n", J++, I.coor_y, I.coor_z);
        I.s_material->print();
    }
}
