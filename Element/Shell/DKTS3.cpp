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

#include "DKTS3.h"
#include <Domain/DomainBase.h>
#include <Element/Plate/DKT3.h>
#include <Material/Material.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>

DKTS3::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const double E, const double F, unique_ptr<Material>&& M)
    : eccentricity(E)
    , factor(F)
    , s_material(std::forward<unique_ptr<Material>>(M)) {}

DKTS3::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const SectionIntegrationPoint& old_obj)
    : eccentricity(old_obj.eccentricity)
    , factor(old_obj.factor)
    , s_material(old_obj.s_material == nullptr ? nullptr : old_obj.s_material->get_copy()) {}

DKTS3::IntegrationPoint::IntegrationPoint(vec&& C)
    : coor(std::forward<vec>(C))
    , BM(3, s_size / 2)
    , BP(3, s_size / 2) {}

mat DKTS3::form_coor(const mat& C) {
    const auto& X1 = C(0, 0),& X2 = C(1, 0),& X3 = C(2, 0),& Y1 = C(0, 1),& Y2 = C(1, 1),& Y3 = C(2, 1);

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

field<mat> DKTS3::form_transform(const mat& C) {
    const auto& X1 = C(0, 0),& X2 = C(1, 0),& X3 = C(2, 0),& Y1 = C(0, 1),& Y2 = C(1, 1),& Y3 = C(2, 1);

    const auto DX4 = X2 - X1, DX5 = X3 - X2, DX6 = X1 - X3;
    const auto DY4 = Y2 - Y1, DY5 = Y3 - Y2, DY6 = Y1 - Y3;

    const auto L4 = sqrt(DX4 * DX4 + DY4 * DY4);
    const auto L5 = sqrt(DX5 * DX5 + DY5 * DY5);
    const auto L6 = sqrt(DX6 * DX6 + DY6 * DY6);

    const auto C4 = DY4 / L4, C5 = DY5 / L5, C6 = DY6 / L6;
    const auto S4 = -DX4 / L4, S5 = -DX5 / L5, S6 = -DX6 / L6;

    mat BMX(6, 9, fill::zeros), BMY(6, 9, fill::zeros);

    BMX(0, 0) = BMX(1, 3) = BMX(2, 6) = 1.;
    BMX(3, 0) = BMX(3, 3) = BMX(4, 3) = BMX(4, 6) = BMX(5, 6) = BMX(5, 0) = .5;
    BMX(3, 2) = -(BMX(3, 5) = .125 * DY4);
    BMX(4, 5) = -(BMX(4, 8) = .125 * DY5);
    BMX(5, 8) = -(BMX(5, 2) = .125 * DY6);

    BMY(0, 1) = BMY(1, 4) = BMY(2, 7) = 1.;
    BMY(3, 1) = BMY(3, 4) = BMY(4, 4) = BMY(4, 7) = BMY(5, 7) = BMY(5, 1) = .5;
    BMY(3, 5) = -(BMY(3, 2) = .125 * DX4);
    BMY(4, 8) = -(BMY(4, 5) = .125 * DX5);
    BMY(5, 2) = -(BMY(5, 8) = .125 * DX6);

    mat BPX(6, 9, fill::zeros), BPY(6, 9, fill::zeros);

    BPX(0, 2) = BPX(1, 5) = BPX(2, 8) = 1.;

    BPX(3, 0) = -(BPX(3, 3) = 1.5 * S4 / L4);
    BPX(3, 4) = BPX(3, 1) = -.75 * C4 * S4;
    BPX(3, 5) = BPX(3, 2) = .5 * C4 * C4 - .25 * S4 * S4;

    BPX(4, 3) = -(BPX(4, 6) = 1.5 * S5 / L5);
    BPX(4, 4) = BPX(4, 7) = -.75 * C5 * S5;
    BPX(4, 5) = BPX(4, 8) = .5 * C5 * C5 - .25 * S5 * S5;

    BPX(5, 6) = -(BPX(5, 0) = 1.5 * S6 / L6);
    BPX(5, 1) = BPX(5, 7) = -.75 * C6 * S6;
    BPX(5, 2) = BPX(5, 8) = .5 * C6 * C6 - .25 * S6 * S6;

    BPY(0, 1) = BPY(1, 4) = BPY(2, 7) = -1.;

    BPY(3, 3) = -(BPY(3, 0) = 1.5 * C4 / L4);
    BPY(3, 1) = BPY(3, 4) = .25 * C4 * C4 - .5 * S4 * S4;
    BPY(3, 2) = BPY(3, 5) = .75 * C4 * S4;

    BPY(4, 6) = -(BPY(4, 3) = 1.5 * C5 / L5);
    BPY(4, 4) = BPY(4, 7) = .25 * C5 * C5 - .5 * S5 * S5;
    BPY(4, 5) = BPY(4, 8) = .75 * C5 * S5;

    BPY(5, 0) = -(BPY(5, 6) = 1.5 * C6 / L6);
    BPY(5, 1) = BPY(5, 7) = .25 * C6 * C6 - .5 * S6 * S6;
    BPY(5, 2) = BPY(5, 8) = .75 * C6 * S6;

    return {BMX, BMY, BPX, BPY};
}

DKTS3::DKTS3(const unsigned T, uvec&& N, const unsigned M, const double TH, const unsigned IP)
    : ShellBase(T, s_node, s_dof, std::forward<uvec>(N), uvec{M}, false, {DOF::U1, DOF::U2, DOF::U3, DOF::UR1, DOF::UR2, DOF::UR3})
    , thickness(TH)
    , num_ip(IP > 20 ? 20 : IP) {}

int DKTS3::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get<Material>(material_tag(0));

    auto& mat_stiff = mat_proto->get_initial_stiffness();

    direction_cosine(); // compute transformation matrix between global and local

    const auto coor = get_local_coordinate();
    const auto ele_coor = form_coor(coor);
    const auto dkt_trans = form_transform(coor);

    const auto area = .5 * det(ele_coor(span(0, 2), span(0, 2)));

    const mat inv_coor = inv(ele_coor);

    const auto& BMX = dkt_trans(0),& BMY = dkt_trans(1),& BPX = dkt_trans(2),& BPY = dkt_trans(3);

    // along thickness
    const IntegrationPlan t_plan(1, num_ip, IntegrationType::GAUSS);

    mat m_stiffness(s_size / 2, s_size / 2, fill::zeros);
    mat p_stiffness(s_size / 2, s_size / 2, fill::zeros);

    int_pt.clear();
    int_pt.reserve(3);
    for(unsigned I = 0; I < 3; ++I) {
        int_pt.emplace_back(vec{ele_coor(I + 3llu, 1), ele_coor(I + 3llu, 2)});

        auto& m_ip = int_pt.back();

        const mat pn_pxy = shape::triangle(m_ip.coor, 1) * inv_coor;
        m_ip.BM.row(0) = pn_pxy.row(0) * BMX;
        m_ip.BM.row(1) = pn_pxy.row(1) * BMY;
        m_ip.BM.row(2) = pn_pxy.row(0) * BMY + pn_pxy.row(1) * BMX;
        m_ip.BP.row(0) = pn_pxy.row(0) * BPX;
        m_ip.BP.row(1) = pn_pxy.row(1) * BPY;
        m_ip.BP.row(2) = pn_pxy.row(0) * BPY + pn_pxy.row(1) * BPX;

        auto& s_ip = m_ip.sec_int_pt;
        s_ip.clear();
        s_ip.reserve(t_plan.n_rows);
        for(unsigned J = 0; J < t_plan.n_rows; ++J) {
            const auto t_eccentricity = .5 * t_plan(J, 0) * thickness;
            s_ip.emplace_back(t_eccentricity, thickness * t_plan(J, 1) * area / 6., mat_proto->get_copy());
            m_stiffness += s_ip.back().factor * m_ip.BM.t() * mat_stiff * m_ip.BM;
            p_stiffness += t_eccentricity * t_eccentricity * s_ip.back().factor * m_ip.BP.t() * mat_stiff * m_ip.BP;
        }
    }

    trial_stiffness = current_stiffness = initial_stiffness = transform_from_local_to_global(reshuffle(m_stiffness, p_stiffness));

    return SUANPAN_SUCCESS;
}

int DKTS3::update_status() {
    // separate displacement vector
    const auto t_disp = transform_from_global_to_local(get_trial_displacement());
    vec m_disp(9), p_disp(9);
    for(unsigned I = 0, J = 0; I < s_size; I += s_dof, J += 3) {
        const span t_span(J, J + 2llu);
        m_disp(t_span) = t_disp(I + m_dof);
        p_disp(t_span) = t_disp(I + p_dof);
    }

    mat t_stiffness(3, 3), m_stiffness(9, 9, fill::zeros), p_stiffness(9, 9, fill::zeros);
    vec t_stress(3), m_resistance(9, fill::zeros), p_resistance(9, fill::zeros);

    for(const auto& I : int_pt) {
        const vec m_strain = I.BM * m_disp, p_strain = I.BP * p_disp;

        for(const auto& J : I.sec_int_pt) if(J.s_material->update_trial_status(m_strain + J.eccentricity * p_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        t_stress.zeros();
        t_stiffness.zeros();
        for(const auto& J : I.sec_int_pt) {
            t_stress += J.factor * J.s_material->get_trial_stress();
            t_stiffness += J.factor * J.s_material->get_trial_stiffness();
        }
        m_resistance += I.BM.t() * t_stress;
        m_stiffness += I.BM.t() * t_stiffness * I.BM;

        t_stress.zeros();
        t_stiffness.zeros();
        for(const auto& J : I.sec_int_pt) {
            t_stress += J.factor * J.eccentricity * J.s_material->get_trial_stress();
            t_stiffness += J.factor * J.eccentricity * J.eccentricity * J.s_material->get_trial_stiffness();
        }
        p_resistance += I.BP.t() * t_stress;
        p_stiffness += I.BP.t() * t_stiffness * I.BP;
    }

    trial_resistance = transform_from_local_to_global(reshuffle(m_resistance, p_resistance));
    trial_stiffness = transform_from_local_to_global(reshuffle(m_stiffness, p_stiffness));

    return SUANPAN_SUCCESS;
}

int DKTS3::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.s_material->commit_status();
    return code;
}

int DKTS3::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.s_material->clear_status();
    return code;
}

int DKTS3::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.s_material->reset_status();
    return code;
}

vector<vec> DKTS3::record(const OutputType P) {
    vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) for(const auto& K : J.s_material->record(P)) data.emplace_back(K);
    return data;
}

void DKTS3::print() { suanpan_info("A three node planar shell element using CST for membrane action and DKT3 for plate action.\n"); }

#ifdef SUANPAN_VTK
#include <vtkTriangle.h>

void DKTS3::Setup() {
    vtk_cell = vtkSmartPointer<vtkTriangle>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < s_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

void DKTS3::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, s_node, fill::zeros);

    if(OutputType::A == type) t_disp = reshape(get_current_acceleration(), s_dof, s_node);
    else if(OutputType::V == type) t_disp = reshape(get_current_velocity(), s_dof, s_node);
    else if(OutputType::U == type) t_disp = reshape(get_current_displacement(), s_dof, s_node);

    for(unsigned I = 0; I < s_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void DKTS3::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat t_mat = amplifier * mat(reshape(get_current_displacement(), s_dof, s_node)).rows(0, 2).t() + get_coordinate(3);
    for(unsigned I = 0; I < s_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), t_mat(I, 0), t_mat(I, 1), t_mat(I, 2));
}

#endif
