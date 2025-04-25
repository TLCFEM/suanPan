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

#include "DKT3.h"
#include <Domain/DomainBase.h>
#include <Material/Material.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>

DKT3::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const double E, const double F, unique_ptr<Material>&& M)
    : eccentricity(E)
    , factor(F)
    , p_material(std::move(M)) {}

DKT3::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const SectionIntegrationPoint& old_obj)
    : eccentricity(old_obj.eccentricity)
    , factor(old_obj.factor)
    , p_material(old_obj.p_material->get_copy()) {}

DKT3::IntegrationPoint::IntegrationPoint(vec&& C)
    : coor(std::move(C))
    , strain_mat(3, p_size) {}

mat DKT3::form_coor(const mat& C) {
    const auto &X1 = C(0, 0), &X2 = C(1, 0), &X3 = C(2, 0);
    const auto &Y1 = C(0, 1), &Y2 = C(1, 1), &Y3 = C(2, 1);

    mat coor(6, 6);

    coor.col(0).fill(1.);
    coor(0, 1) = X1;
    coor(1, 1) = X2;
    coor(2, 1) = X3;
    coor(3, 1) = .5 * (X1 + X2);
    coor(4, 1) = .5 * (X2 + X3);
    coor(5, 1) = .5 * (X3 + X1);
    coor(0, 2) = Y1;
    coor(1, 2) = Y2;
    coor(2, 2) = Y3;
    coor(3, 2) = .5 * (Y1 + Y2);
    coor(4, 2) = .5 * (Y2 + Y3);
    coor(5, 2) = .5 * (Y3 + Y1);
    coor.col(3) = coor.col(1) % coor.col(2);
    coor.col(4) = square(coor.col(1));
    coor.col(5) = square(coor.col(2));

    return coor;
}

field<mat> DKT3::form_transform(const mat& C) {
    const auto &X1 = C(0, 0), &X2 = C(1, 0), &X3 = C(2, 0), &Y1 = C(0, 1), &Y2 = C(1, 1), &Y3 = C(2, 1);

    const auto DX4 = X2 - X1, DX5 = X3 - X2, DX6 = X1 - X3;
    const auto DY4 = Y2 - Y1, DY5 = Y3 - Y2, DY6 = Y1 - Y3;

    const auto L4 = sqrt(DX4 * DX4 + DY4 * DY4);
    const auto L5 = sqrt(DX5 * DX5 + DY5 * DY5);
    const auto L6 = sqrt(DX6 * DX6 + DY6 * DY6);

    const auto C4 = DY4 / L4, C5 = DY5 / L5, C6 = DY6 / L6;
    const auto S4 = -DX4 / L4, S5 = -DX5 / L5, S6 = -DX6 / L6;

    mat BX(6, p_size, fill::zeros), BY(6, p_size, fill::zeros);

    BX(0, 2) = BX(1, 5) = BX(2, 8) = 1.;

    BX(3, 0) = -(BX(3, 3) = 1.5 * S4 / L4);
    BX(3, 4) = BX(3, 1) = -.75 * C4 * S4;
    BX(3, 5) = BX(3, 2) = .5 * C4 * C4 - .25 * S4 * S4;

    BX(4, 3) = -(BX(4, 6) = 1.5 * S5 / L5);
    BX(4, 4) = BX(4, 7) = -.75 * C5 * S5;
    BX(4, 5) = BX(4, 8) = .5 * C5 * C5 - .25 * S5 * S5;

    BX(5, 6) = -(BX(5, 0) = 1.5 * S6 / L6);
    BX(5, 1) = BX(5, 7) = -.75 * C6 * S6;
    BX(5, 2) = BX(5, 8) = .5 * C6 * C6 - .25 * S6 * S6;

    BY(0, 1) = BY(1, 4) = BY(2, 7) = -1.;

    BY(3, 3) = -(BY(3, 0) = 1.5 * C4 / L4);
    BY(3, 1) = BY(3, 4) = .25 * C4 * C4 - .5 * S4 * S4;
    BY(3, 2) = BY(3, 5) = .75 * C4 * S4;

    BY(4, 6) = -(BY(4, 3) = 1.5 * C5 / L5);
    BY(4, 4) = BY(4, 7) = .25 * C5 * C5 - .5 * S5 * S5;
    BY(4, 5) = BY(4, 8) = .75 * C5 * S5;

    BY(5, 0) = -(BY(5, 6) = 1.5 * C6 / L6);
    BY(5, 1) = BY(5, 7) = .25 * C6 * C6 - .5 * S6 * S6;
    BY(5, 2) = BY(5, 8) = .75 * C6 * S6;

    return {BX, BY};
}

DKT3::DKT3(const unsigned T, uvec&& NT, const unsigned MT, const double TH, const unsigned IPN)
    : MaterialElement2D(T, p_node, p_dof, std::move(NT), uvec{MT}, false, {DOF::U1, DOF::U2, DOF::UR3})
    , thickness(TH)
    , num_section_ip(IPN) {}

int DKT3::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get<Material>(material_tag(0));

    auto& ini_stiffness = mat_proto->get_initial_stiffness();

    const auto coor = get_coordinate(2);
    const auto ele_coor = form_coor(coor);
    const auto trans_mat = form_transform(coor);

    const auto area = .5 * det(ele_coor(span(0, 2), span(0, 2)));

    const mat inv_coor = inv(ele_coor);

    const auto& BX = trans_mat(0);
    const auto& BY = trans_mat(1);

    const IntegrationPlan sec_plan(1, num_section_ip, IntegrationType::GAUSS);

    initial_stiffness.zeros(p_size, p_size);

    int_pt.clear();
    int_pt.reserve(3);
    for(auto I = 0; I < 3; ++I) {
        int_pt.emplace_back(vec{ele_coor(I + 3llu, 1), ele_coor(I + 3llu, 2)});

        auto& c_pt = int_pt.back();
        auto& strain_mat = c_pt.strain_mat;

        const mat pn_pxy = shape::triangle(c_pt.coor, 1) * inv_coor;
        strain_mat.row(0) = pn_pxy.row(0) * BX;
        strain_mat.row(1) = pn_pxy.row(1) * BY;
        strain_mat.row(2) = pn_pxy.row(0) * BY + pn_pxy.row(1) * BX;

        auto& c_ip = c_pt.sec_int_pt;
        c_ip.clear();
        c_ip.reserve(num_section_ip);
        for(unsigned J = 0; J < num_section_ip; ++J) {
            const auto t_eccentricity = .5 * sec_plan(J, 0) * thickness;
            c_ip.emplace_back(t_eccentricity, thickness * sec_plan(J, 1) * area / 6., mat_proto->get_copy());
            initial_stiffness += t_eccentricity * t_eccentricity * c_ip.back().factor * strain_mat.t() * ini_stiffness * strain_mat;
        }
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int DKT3::update_status() {
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

int DKT3::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.p_material->clear_status();
    return code;
}

int DKT3::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.p_material->commit_status();
    return code;
}

int DKT3::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.p_material->reset_status();
    return code;
}

std::vector<vec> DKT3::record(const OutputType P) {
    std::vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) append_to(data, J.p_material->record(P));
    return data;
}

void DKT3::print() {
    suanpan_info("A DKT triangular plate element connects:", node_encoding);
}

#ifdef SUANPAN_VTK
#include <vtkTriangle.h>

void DKT3::Setup() {
    vtk_cell = vtkSmartPointer<vtkTriangle>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < p_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void DKT3::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, p_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(2, 4) = reshape(get_current_acceleration(), p_dof, p_node);
    else if(OutputType::V == type) t_disp.rows(2, 4) = reshape(get_current_velocity(), p_dof, p_node);
    else if(OutputType::U == type) t_disp.rows(2, 4) = reshape(get_current_displacement(), p_dof, p_node);

    for(unsigned I = 0; I < p_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void DKT3::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const auto ele_coor = get_coordinate(2);
    const mat ele_disp = reshape(get_current_displacement(), p_dof, p_node);
    for(unsigned I = 0; I < p_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_coor(I, 0), ele_coor(I, 1), amplifier * ele_disp(0, I));
}

#endif
