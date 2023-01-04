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

#include "DKT4.h"
#include <Domain/DomainBase.h>
#include <Material/Material.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>

DKT4::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const double E, const double F, unique_ptr<Material>&& M)
    : eccentricity(E)
    , factor(F)
    , p_material(std::forward<unique_ptr<Material>>(M)) {}

DKT4::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const SectionIntegrationPoint& old_obj)
    : eccentricity(old_obj.eccentricity)
    , factor(old_obj.factor)
    , p_material(old_obj.p_material->get_copy()) {}

DKT4::IntegrationPoint::IntegrationPoint(vec&& C)
    : coor(std::forward<vec>(C))
    , strain_mat(3, p_size) {}

field<mat> DKT4::form_transform(const mat& C) {
    const auto &X1 = C(0, 0), &X2 = C(1, 0), &X3 = C(2, 0), &X4 = C(3, 0);
    const auto &Y1 = C(0, 1), &Y2 = C(1, 1), &Y3 = C(2, 1), &Y4 = C(3, 1);

    const auto DX5 = X2 - X1, DX6 = X3 - X2, DX7 = X4 - X3, DX8 = X1 - X4;
    const auto DY5 = Y2 - Y1, DY6 = Y3 - Y2, DY7 = Y4 - Y3, DY8 = Y1 - Y4;

    const auto L5 = sqrt(DX5 * DX5 + DY5 * DY5);
    const auto L6 = sqrt(DX6 * DX6 + DY6 * DY6);
    const auto L7 = sqrt(DX7 * DX7 + DY7 * DY7);
    const auto L8 = sqrt(DX8 * DX8 + DY8 * DY8);

    const auto C5 = DY5 / L5, C6 = DY6 / L6, C7 = DY7 / L7, C8 = DY8 / L8;
    const auto S5 = -DX5 / L5, S6 = -DX6 / L6, S7 = -DX7 / L7, S8 = -DX8 / L8;

    mat BX(8, p_size, fill::zeros), BY(8, p_size, fill::zeros);

    BX(0, 2) = BX(1, 5) = BX(2, 8) = BX(3, 11) = 1.;

    BX(4, 0) = -(BX(4, 3) = 1.5 * S5 / L5);
    BX(4, 4) = BX(4, 1) = -.75 * C5 * S5;
    BX(4, 5) = BX(4, 2) = .5 * C5 * C5 - .25 * S5 * S5;

    BX(5, 3) = -(BX(5, 6) = 1.5 * S6 / L6);
    BX(5, 7) = BX(5, 4) = -.75 * C6 * S6;
    BX(5, 8) = BX(5, 5) = .5 * C6 * C6 - .25 * S6 * S6;

    BX(6, 6) = -(BX(6, 9) = 1.5 * S7 / L7);
    BX(6, 10) = BX(6, 7) = -.75 * C7 * S7;
    BX(6, 11) = BX(6, 8) = .5 * C7 * C7 - .25 * S7 * S7;

    BX(7, 9) = -(BX(7, 0) = 1.5 * S8 / L8);
    BX(7, 1) = BX(7, 10) = -.75 * C8 * S8;
    BX(7, 2) = BX(7, 11) = .5 * C8 * C8 - .25 * S8 * S8;

    BY(0, 1) = BY(1, 4) = BY(2, 7) = BY(3, 10) = -1.;

    BY(4, 3) = -(BY(4, 0) = 1.5 * C5 / L5);
    BY(4, 4) = BY(4, 1) = .25 * C5 * C5 - .5 * S5 * S5;
    BY(4, 5) = BY(4, 2) = .75 * C5 * S5;

    BY(5, 6) = -(BY(5, 3) = 1.5 * C6 / L6);
    BY(5, 7) = BY(5, 4) = .25 * C6 * C6 - .5 * S6 * S6;
    BY(5, 8) = BY(5, 5) = .75 * C6 * S6;

    BY(6, 9) = -(BY(6, 6) = 1.5 * C7 / L7);
    BY(6, 10) = BY(6, 7) = .25 * C7 * C7 - .5 * S7 * S7;
    BY(6, 11) = BY(6, 8) = .75 * C7 * S7;

    BY(7, 0) = -(BY(7, 9) = 1.5 * C8 / L8);
    BY(7, 1) = BY(7, 10) = .25 * C8 * C8 - .5 * S8 * S8;
    BY(7, 2) = BY(7, 11) = .75 * C8 * S8;

    return {BX, BY};
}

DKT4::DKT4(const unsigned T, uvec&& NT, const unsigned MT, const double TH, const unsigned IPN)
    : MaterialElement2D(T, p_node, p_dof, std::forward<uvec>(NT), uvec{MT}, false, {DOF::U1, DOF::U2, DOF::UR3})
    , thickness(TH)
    , num_section_ip(IPN) {}

int DKT4::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get<Material>(material_tag(0));

    auto& ini_stiffness = mat_proto->get_initial_stiffness();

    const auto coor = get_coordinate(2);
    const auto trans_mat = form_transform(coor);

    const auto& BX = trans_mat(0);
    const auto& BY = trans_mat(1);

    mat ele_coor(8, 2);
    ele_coor.rows(0, 3) = coor;
    ele_coor.row(4) = .5 * (coor.row(0) + coor.row(1));
    ele_coor.row(5) = .5 * (coor.row(1) + coor.row(2));
    ele_coor.row(6) = .5 * (coor.row(2) + coor.row(3));
    ele_coor.row(7) = .5 * (coor.row(3) + coor.row(0));

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);
    const IntegrationPlan sec_plan(1, num_section_ip, IntegrationType::GAUSS);

    initial_stiffness.zeros(p_size, p_size);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        int_pt.emplace_back(vec{plan(I, 0), plan(I, 1)});

        auto& c_pt = int_pt.back();

        const auto pn = shape::quad(c_pt.coor, 1, 8);
        const mat jacob = pn * ele_coor;
        const mat pn_pxy = solve(jacob, pn);

        auto& strain_mat = c_pt.strain_mat;
        strain_mat.row(0) = pn_pxy.row(0) * BX;
        strain_mat.row(1) = pn_pxy.row(1) * BY;
        strain_mat.row(2) = pn_pxy.row(0) * BY + pn_pxy.row(1) * BX;

        const auto t_weight = .5 * thickness * plan(I, 2) * det(jacob);
        auto& c_ip = c_pt.sec_int_pt;
        c_ip.clear();
        c_ip.reserve(num_section_ip);
        for(unsigned J = 0; J < num_section_ip; ++J) {
            const auto t_eccentricity = .5 * sec_plan(J, 0) * thickness;
            c_ip.emplace_back(t_eccentricity, t_weight * sec_plan(J, 1), mat_proto->get_copy());
            initial_stiffness += t_eccentricity * t_eccentricity * c_ip.back().factor * strain_mat.t() * ini_stiffness * strain_mat;
        }
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int DKT4::update_status() {
    const auto trial_disp = get_trial_displacement();

    trial_resistance.zeros(p_size);
    trial_stiffness.zeros(p_size, p_size);
    for(const auto& I : int_pt) {
        const vec p_strain = I.strain_mat * trial_disp;
        for(const auto& J : I.sec_int_pt) {
            if(J.p_material->update_trial_status(J.eccentricity * p_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
            trial_stiffness += J.eccentricity * J.eccentricity * J.factor * I.strain_mat.t() * J.p_material->get_trial_stiffness() * I.strain_mat;
            trial_resistance += J.eccentricity * J.factor * I.strain_mat.t() * J.p_material->get_trial_stress();
        }
    }

    return SUANPAN_SUCCESS;
}

int DKT4::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.p_material->clear_status();
    return code;
}

int DKT4::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.p_material->commit_status();
    return code;
}

int DKT4::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.p_material->reset_status();
    return code;
}

vector<vec> DKT4::record(const OutputType P) {
    vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) for(const auto& K : J.p_material->record(P)) data.emplace_back(K);
    return data;
}

void DKT4::print() {
    suanpan_info("A DKT quadrilateral plate element connects:", node_encoding);
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void DKT4::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < p_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void DKT4::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, p_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(2, 4) = reshape(get_current_acceleration(), p_dof, p_node);
    else if(OutputType::V == type) t_disp.rows(2, 4) = reshape(get_current_velocity(), p_dof, p_node);
    else if(OutputType::U == type) t_disp.rows(2, 4) = reshape(get_current_displacement(), p_dof, p_node);

    for(unsigned I = 0; I < p_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void DKT4::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const auto ele_coor = get_coordinate(2);
    const mat ele_disp = reshape(get_current_displacement(), p_dof, p_node);
    for(unsigned I = 0; I < p_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_coor(I, 0), ele_coor(I, 1), amplifier * ele_disp(0, I));
}

#endif
