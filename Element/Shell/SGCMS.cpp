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

#include "SGCMS.h"

#include <Domain/DomainBase.h>
#include <Material/Material.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>

const mat SGCMS::mapping = [] {
    mat t_mapping(4, 4);
    t_mapping.fill(.25);
    t_mapping(1, 0) = t_mapping(1, 3) = t_mapping(2, 0) = t_mapping(2, 1) = t_mapping(3, 1) = t_mapping(3, 3) = -.25;
    return t_mapping;
}();

SGCMS::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const double E, const double F, unique_ptr<Material>&& M)
    : eccentricity(E)
    , factor(F)
    , s_material(std::move(M)) {}

SGCMS::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const SectionIntegrationPoint& old_obj)
    : eccentricity(old_obj.eccentricity)
    , factor(old_obj.factor)
    , s_material(nullptr == old_obj.s_material ? nullptr : old_obj.s_material->get_copy()) {}

SGCMS::IntegrationPoint::IntegrationPoint(vec&& C)
    : coor(std::move(C))
    , BP(3, 12, fill::zeros) {}

field<mat> SGCMS::form_plate_transformation(const mat& C) {
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

    mat BX(8, 12, fill::zeros), BY(8, 12, fill::zeros);

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

mat SGCMS::form_drilling_n(const vec& coor, const vec& lxy) {
    mat poly(2, 12, fill::zeros);

    auto &X = coor(0), &Y = coor(1);

    auto &LX1 = lxy(0), &LX2 = lxy(1), &LX3 = lxy(2), &LX4 = lxy(3);
    auto &LY1 = lxy(4), &LY2 = lxy(5), &LY3 = lxy(6), &LY4 = lxy(7);

    const auto XX = X * X - 1., YY = Y * Y - 1., YP = 1. + Y, YM = 1. - Y, XP = 1. + X, XM = 1. - X;

    poly(0, 2) = LX1 * XX * YM - LX4 * YY * XM;
    poly(0, 5) = LX2 * YY * XP - LX1 * XX * YM;
    poly(0, 8) = LX3 * XX * YP - LX2 * YY * XP;
    poly(0, 11) = LX4 * YY * XM - LX3 * XX * YP;
    poly(1, 2) = LY1 * XX * YM - LY4 * YY * XM;
    poly(1, 5) = LY2 * YY * XP - LY1 * XX * YM;
    poly(1, 8) = LY3 * XX * YP - LY2 * YY * XP;
    poly(1, 11) = LY4 * YY * XM - LY3 * XX * YP;

    return poly /= 16.;
}

mat SGCMS::form_drilling_dn(const vec& coor, const vec& lxy) {
    mat poly(2, 8);

    auto &X = coor(0), &Y = coor(1);

    auto &LX1 = lxy(0), &LX2 = lxy(1), &LX3 = lxy(2), &LX4 = lxy(3);
    auto &LY1 = lxy(4), &LY2 = lxy(5), &LY3 = lxy(6), &LY4 = lxy(7);

    const auto X2 = 2. * X, Y2 = 2. * Y, XP = X + 1., XM = X - 1., YP = Y + 1., YM = Y - 1.;

    poly(0, 0) = +YM * (LX4 * YP - LX1 * X2);
    poly(0, 1) = +YM * (LX2 * YP + LX1 * X2);
    poly(0, 2) = -YP * (LX2 * YM - LX3 * X2);
    poly(0, 3) = -YP * (LX4 * YM + LX3 * X2);
    poly(0, 4) = +YM * (LY4 * YP - LY1 * X2);
    poly(0, 5) = +YM * (LY2 * YP + LY1 * X2);
    poly(0, 6) = -YP * (LY2 * YM - LY3 * X2);
    poly(0, 7) = -YP * (LY4 * YM + LY3 * X2);
    poly(1, 0) = -XM * (LX1 * XP - LX4 * Y2);
    poly(1, 1) = +XP * (LX1 * XM + LX2 * Y2);
    poly(1, 2) = +XP * (LX3 * XM - LX2 * Y2);
    poly(1, 3) = -XM * (LX3 * XP + LX4 * Y2);
    poly(1, 4) = -XM * (LY1 * XP - LY4 * Y2);
    poly(1, 5) = +XP * (LY1 * XM + LY2 * Y2);
    poly(1, 6) = +XP * (LY3 * XM - LY2 * Y2);
    poly(1, 7) = -XM * (LY3 * XP + LY4 * Y2);

    return poly /= 16.;
}

mat SGCMS::form_displacement_dn(const mat& pn_pxy, const mat& pnt_pxy) {
    mat poly(3, 12, fill::zeros);

    for(unsigned J = 0, K = 0, L = 1, M = 2, N = 4; J < s_node; ++J, K += 3, L += 3, M += 3, ++N) {
        poly(0, K) = poly(2, L) = pn_pxy(0, J);
        poly(2, K) = poly(1, L) = pn_pxy(1, J);
        poly(0, M) = pnt_pxy(0, J);
        poly(1, M) = pnt_pxy(1, N);
        poly(2, M) = pnt_pxy(0, N) + pnt_pxy(1, J);
    }

    return poly;
}

SGCMS::SGCMS(const unsigned T, uvec&& N, const unsigned M, const double TH, const bool NL)
    : ShellBase(T, s_node, s_dof, std::move(N), uvec{M}, NL, {DOF::U1, DOF::U2, DOF::U3, DOF::UR1, DOF::UR2, DOF::UR3})
    , thickness(TH) {}

int SGCMS::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get<Material>(material_tag(0));

    direction_cosine();

    auto& mat_stiff = mat_proto->get_initial_stiffness();

    auto ele_coor = get_local_coordinate();

    // in-plane
    const IntegrationPlan m_plan(2, 2, IntegrationType::IRONS);

    //
    // membrane action
    //
    vec diff_coor(8);
    diff_coor(0) = ele_coor(1, 1) - ele_coor(0, 1);
    diff_coor(1) = ele_coor(2, 1) - ele_coor(1, 1);
    diff_coor(2) = ele_coor(3, 1) - ele_coor(2, 1);
    diff_coor(3) = ele_coor(0, 1) - ele_coor(3, 1);
    diff_coor(4) = ele_coor(0, 0) - ele_coor(1, 0);
    diff_coor(5) = ele_coor(1, 0) - ele_coor(2, 0);
    diff_coor(6) = ele_coor(2, 0) - ele_coor(3, 0);
    diff_coor(7) = ele_coor(3, 0) - ele_coor(0, 0);

    const mat iso_mapping = trans(mapping * ele_coor);

    mat N(11, 12, fill::zeros), H(11, 11, fill::zeros), HT(11, 11, fill::zeros);

    int_pt.clear();
    int_pt.reserve(m_plan.n_rows);
    for(unsigned I = 0; I < m_plan.n_rows; ++I) {
        const auto &X = m_plan(I, 0), &Y = m_plan(I, 1);

        int_pt.emplace_back(vec{X, Y});

        auto& c_ip = int_pt.back();

        const auto pn = shape::quad(c_ip.coor, 1);
        const mat jacob = pn * ele_coor;

        const auto poly_stress = shape::stress11(iso_mapping * vec{0., X, Y, X * Y});
        c_ip.BM = solve(mat_stiff, poly_stress);

        const auto factor = det(jacob) * m_plan(I, 2) * thickness;
        N += factor * poly_stress.t() * form_displacement_dn(solve(jacob, pn), solve(jacob, form_drilling_dn(c_ip.coor, diff_coor)));
        H += factor * poly_stress.t() * c_ip.BM;
        HT += factor * c_ip.BM.t() * mat_stiff * c_ip.BM;
    }

    const mat NT = solve(H, N);

    //
    // plate action
    //
    const auto trans = form_plate_transformation(ele_coor);
    const auto &BX = trans(0), &BY = trans(1);

    ele_coor.resize(8, 2);
    ele_coor.row(4) = .5 * (ele_coor.row(0) + ele_coor.row(1));
    ele_coor.row(5) = .5 * (ele_coor.row(1) + ele_coor.row(2));
    ele_coor.row(6) = .5 * (ele_coor.row(2) + ele_coor.row(3));
    ele_coor.row(7) = .5 * (ele_coor.row(3) + ele_coor.row(0));

    // along thickness
    const IntegrationPlan t_plan(1, 3, IntegrationType::GAUSS);

    mat::fixed<12, 12> p_stiffness(fill::zeros), mp_stiffness(fill::zeros), pm_stiffness(fill::zeros);

    for(unsigned I = 0; I < m_plan.n_rows; ++I) {
        auto& c_pt = int_pt.at(I);

        // update membrane strain matrix first
        c_pt.BM *= NT;

        const auto pn = shape::quad(c_pt.coor, 1, 8);
        const mat jacob = pn * ele_coor;
        const mat pn_pxy = solve(jacob, pn);

        c_pt.BP.row(0) = pn_pxy.row(0) * BX;
        c_pt.BP.row(1) = pn_pxy.row(1) * BY;
        c_pt.BP.row(2) = pn_pxy.row(0) * BY + pn_pxy.row(1) * BX;

        const auto t_weight = .5 * thickness * m_plan(I, 2) * det(jacob);
        auto& c_ip = c_pt.sec_int_pt;
        c_ip.clear();
        c_ip.reserve(t_plan.n_rows);
        for(unsigned J = 0; J < t_plan.n_rows; ++J) {
            const auto t_eccentricity = .5 * t_plan(J, 0) * thickness;
            c_ip.emplace_back(t_eccentricity, t_weight * t_plan(J, 1), mat_proto->get_copy());
            p_stiffness += c_pt.BP.t() * mat_stiff * c_pt.BP * c_ip.back().factor * pow(t_eccentricity, 2.);
            mp_stiffness += c_pt.BM.t() * mat_stiff * c_pt.BP * c_ip.back().factor * t_eccentricity;
            pm_stiffness += c_pt.BP.t() * mat_stiff * c_pt.BM * c_ip.back().factor * t_eccentricity;
        }
    }

    transform_from_local_to_global(initial_stiffness = reshuffle(NT.t() * HT * NT, p_stiffness, mp_stiffness, pm_stiffness));
    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int SGCMS::update_status() {
    direction_cosine();

    // separate displacement vector
    const auto g_disp = get_trial_displacement();
    const auto t_disp = transform_from_global_to_local(g_disp);
    vec::fixed<12> m_disp, p_disp;
    for(unsigned I = 0, J = 0; I < s_size; I += s_dof, J += 3) {
        const span t_span(J, J + 2llu);
        m_disp(t_span) = t_disp(I + m_dof);
        p_disp(t_span) = t_disp(I + p_dof);
    }

    mat33 t_stiffness;
    mat::fixed<12, 12> m_stiffness(fill::zeros), p_stiffness(fill::zeros), mp_stiffness(fill::zeros), pm_stiffness(fill::zeros);
    vec3 t_stress;
    vec::fixed<12> m_resistance(fill::zeros), p_resistance(fill::zeros);

    for(const auto& I : int_pt) {
        const vec m_strain = I.BM * m_disp, p_strain = I.BP * p_disp;

        for(const auto& J : I.sec_int_pt)
            if(J.s_material->update_trial_status(m_strain + J.eccentricity * p_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

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

        t_stiffness.zeros();
        for(const auto& J : I.sec_int_pt) t_stiffness += J.factor * J.eccentricity * J.s_material->get_trial_stiffness();
        mp_stiffness += I.BM.t() * t_stiffness * I.BP;
        pm_stiffness += I.BP.t() * t_stiffness * I.BM;
    }

    trial_resistance = reshuffle(m_resistance, p_resistance);
    trial_stiffness = reshuffle(m_stiffness, p_stiffness, mp_stiffness, pm_stiffness);

    if(is_nlgeom()) trial_geometry = transform_to_global_geometry(trial_stiffness, trial_resistance, g_disp);

    transform_from_local_to_global(trial_resistance);
    transform_from_local_to_global(trial_stiffness);

    return SUANPAN_SUCCESS;
}

int SGCMS::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt)
        for(const auto& J : I.sec_int_pt) code += J.s_material->commit_status();
    return code;
}

int SGCMS::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt)
        for(const auto& J : I.sec_int_pt) code += J.s_material->clear_status();
    return code;
}

int SGCMS::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt)
        for(const auto& J : I.sec_int_pt) code += J.s_material->reset_status();
    return code;
}

std::vector<vec> SGCMS::record(const OutputType P) {
    std::vector<vec> data;
    for(const auto& I : int_pt)
        for(const auto& J : I.sec_int_pt) append_to(data, J.s_material->record(P));
    return data;
}

void SGCMS::print() {
    suanpan_info("A four-node planar shell element using SGCMQ for membrane action and DKT4 for plate action.\n");
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void SGCMS::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < s_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

void SGCMS::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, s_node, fill::zeros);

    if(OutputType::A == type) t_disp = reshape(get_current_acceleration(), s_dof, s_node);
    else if(OutputType::V == type) t_disp = reshape(get_current_velocity(), s_dof, s_node);
    else if(OutputType::U == type) t_disp = reshape(get_current_displacement(), s_dof, s_node);

    for(unsigned I = 0; I < s_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void SGCMS::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat t_mat = amplifier * mat(reshape(get_current_displacement(), s_dof, s_node)).rows(0, 2).t() + get_coordinate(3);
    for(unsigned I = 0; I < s_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), t_mat(I, 0), t_mat(I, 1), t_mat(I, 2));
}

#endif
