/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
    : Section3D(T, MT, TFW * TFT + BFW * BFT + WH * WT, std::move(EC))
    , top_flange_width(TFW)
    , top_flange_thickness(TFT)
    , bottom_flange_width(BFW)
    , bottom_flange_thickness(BFT)
    , web_height(WH)
    , web_thickness(WT)
    , int_pt_num(IP > 20 ? 20 : IP) {}

ISection3D::ISection3D(const unsigned T, vec&& D, const unsigned MT, const unsigned IP, vec&& EC)
    : Section3D(T, MT, D(0) * D(1) + D(2) * D(3) + D(4) * D(5), std::move(EC))
    , top_flange_width(D(0))
    , top_flange_thickness(D(1))
    , bottom_flange_width(D(2))
    , bottom_flange_thickness(D(3))
    , web_height(D(4))
    , web_thickness(D(5))
    , int_pt_num(IP > 20 ? 20 : IP) {}

int ISection3D::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get_material(material_tag);

    access::rw(linear_density) = mat_proto->get_density() * area;

    const auto web_area = web_height * web_thickness;
    const auto b_flange_area = bottom_flange_width * bottom_flange_thickness;
    const auto t_flange_area = top_flange_width * top_flange_thickness;

    const IntegrationPlan plan_ft(1, std::max(bottom_flange_thickness, top_flange_thickness) > .1 * web_height ? 2 : 1, IntegrationPlan::Type::GAUSS);
    const IntegrationPlan plan_wt(1, web_thickness > .1 * std::max(bottom_flange_width, top_flange_width) ? 2 : 1, IntegrationPlan::Type::GAUSS);

    const IntegrationPlan plan_f(1, int_pt_num, IntegrationPlan::Type::GAUSS);
    const auto& plan_w = plan_f;

    int_pt.clear();
    int_pt.reserve(plan_w.n_rows * plan_wt.n_rows + 2u * plan_f.n_rows * plan_ft.n_rows);
    for(auto I = 0u; I < plan_w.n_rows; ++I)
        for(auto J = 0u; J < plan_wt.n_rows; ++J) int_pt.emplace_back(.5 * plan_w(I, 0) * web_height, .5 * plan_wt(J, 0) * web_thickness, .25 * plan_w(I, 1) * plan_wt(J, 1) * web_area, mat_proto->unique_copy());
    if(b_flange_area != 0.)
        for(auto I = 0u; I < plan_f.n_rows; ++I)
            for(auto J = 0u; J < plan_ft.n_rows; ++J) int_pt.emplace_back((.5 * plan_ft(J, 0) + .5) * bottom_flange_thickness + .5 * web_height, .5 * plan_f(I, 0) * bottom_flange_width, .25 * plan_f(I, 1) * plan_ft(J, 1) * b_flange_area, mat_proto->unique_copy());
    if(t_flange_area != 0.)
        for(auto I = 0u; I < plan_f.n_rows; ++I)
            for(auto J = 0u; J < plan_ft.n_rows; ++J) int_pt.emplace_back((-.5 * plan_ft(J, 0) - .5) * top_flange_thickness - .5 * web_height, .5 * plan_f(I, 0) * top_flange_width, .25 * plan_f(I, 1) * plan_ft(J, 1) * t_flange_area, mat_proto->unique_copy());

    initialize_stiffness();

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> ISection3D::unique_copy() { return std::make_unique<ISection3D>(*this); }

void ISection3D::print() {
    suanpan_info("A 3D I-shape section with following integration points.\n");
    auto J = 1;
    for(const auto& I : int_pt) {
        suanpan_info("IP {}: {:.4E}, {:.4E}.\n", J++, I.coor_y, I.coor_z);
        I.s_material->print();
    }
}
