/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "S4.h"
#include <Domain/DomainBase.h>
#include <Material/Material.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>
#include <Toolbox/utility.h>

S4::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const double E, const double F, unique_ptr<Material>&& M)
    : eccentricity(E)
    , factor(F)
    , s_material(std::move(M)) {}

S4::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const SectionIntegrationPoint& old_obj)
    : eccentricity(old_obj.eccentricity)
    , factor(old_obj.factor)
    , s_material(old_obj.s_material == nullptr ? nullptr : old_obj.s_material->get_copy()) {}

S4::IntegrationPoint::IntegrationPoint(vec&& C)
    : coor(std::move(C))
    , BM(3, s_size / 2, fill::zeros)
    , BP(3, s_size / 2, fill::zeros) {}

S4::S4(const unsigned T, uvec&& N, const unsigned M, const double TH, const bool NL)
    : ShellBase(T, s_node, s_dof, std::move(N), uvec{M}, NL, {DOF::U1, DOF::U2, DOF::U3, DOF::UR1, DOF::UR2, DOF::UR3})
    , thickness(TH) {}

int S4::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get<Material>(material_tag(0));

    direction_cosine();

    const auto ele_coor = get_local_coordinate();

    const auto LX1 = ele_coor(1, 1) - ele_coor(0, 1);
    const auto LX2 = ele_coor(2, 1) - ele_coor(1, 1);
    const auto LX3 = ele_coor(3, 1) - ele_coor(2, 1);
    const auto LX4 = ele_coor(0, 1) - ele_coor(3, 1);
    const auto LY1 = ele_coor(0, 0) - ele_coor(1, 0);
    const auto LY2 = ele_coor(1, 0) - ele_coor(2, 0);
    const auto LY3 = ele_coor(2, 0) - ele_coor(3, 0);
    const auto LY4 = ele_coor(3, 0) - ele_coor(0, 0);

    auto& mat_stiff = mat_proto->get_initial_stiffness();

    // Mindlin plate
    // check if proper shear modulus is available
    // not very vital as for multiplier any large value can be chosen
    const auto shear_modulus = mat_proto->get_parameter(ParameterType::SHEARMODULUS);
    if(suanpan::approx_equal(shear_modulus, 0.)) {
        suanpan_error("A zero shear modulus is detected.\n");
        return SUANPAN_FAIL;
    };

    // reduced integration for the Kirchhoff constraint
    vec t_vec(2, fill::zeros);
    const auto n = shape::quad(t_vec, 0);
    auto pn = shape::quad(t_vec, 1);
    mat jacob = pn * ele_coor, pn_pxy = solve(jacob, pn);
    mat penalty_mat(2, 3llu * s_node, fill::zeros);
    for(unsigned I = 0, J = 0; I < 4; ++I, J += 3) {
        penalty_mat(0, J) = pn_pxy(1, I);
        penalty_mat(1, J) = pn_pxy(0, I);
        penalty_mat(0, J + 1llu) = -(penalty_mat(1, J + 2llu) = n(I));
    }
    auto p_stiffness = penalty_stiffness = 10. / 3. * shear_modulus * thickness * det(jacob) * penalty_mat.t() * penalty_mat;

    // in-plane
    const IntegrationPlan m_plan(2, 2, IntegrationType::GAUSS);
    // along thickness
    const IntegrationPlan t_plan(1, 3, IntegrationType::GAUSS);

    mat m_stiffness(s_size / 2, s_size / 2, fill::zeros), pnt(2, 8);

    int_pt.clear();
    int_pt.reserve(m_plan.n_rows);
    for(unsigned I = 0; I < m_plan.n_rows; ++I) {
        int_pt.emplace_back(vec{m_plan(I, 0), m_plan(I, 1)});

        auto& m_ip = int_pt.back();

        pn = shape::quad(m_ip.coor, 1);
        jacob = pn * ele_coor;
        pn_pxy = solve(jacob, pn);
        const auto det_jacob = det(jacob);

        const auto X = 2. * m_plan(I, 0);
        const auto Y = 2. * m_plan(I, 1);

        const auto AA = m_plan(I, 0) + 1.;
        const auto BB = m_plan(I, 0) - 1.;
        const auto CC = m_plan(I, 1) + 1.;
        const auto DD = m_plan(I, 1) - 1.;

        pnt(0, 0) = DD * (LX4 * CC - LX1 * X);
        pnt(0, 1) = DD * (LX2 * CC + LX1 * X);
        pnt(0, 2) = CC * (LX3 * X - LX2 * DD);
        pnt(0, 3) = -CC * (LX3 * X + LX4 * DD);
        pnt(0, 4) = DD * (LY4 * CC - LY1 * X);
        pnt(0, 5) = DD * (LY2 * CC + LY1 * X);
        pnt(0, 6) = CC * (LY3 * X - LY2 * DD);
        pnt(0, 7) = -CC * (LY3 * X + LY4 * DD);
        pnt(1, 0) = BB * (LX4 * Y - LX1 * AA);
        pnt(1, 1) = AA * (LX1 * BB + LX2 * Y);
        pnt(1, 2) = AA * (LX3 * BB - LX2 * Y);
        pnt(1, 3) = -BB * (LX3 * AA + LX4 * Y);
        pnt(1, 4) = BB * (LY4 * Y - LY1 * AA);
        pnt(1, 5) = AA * (LY1 * BB + LY2 * Y);
        pnt(1, 6) = AA * (LY3 * BB - LY2 * Y);
        pnt(1, 7) = -BB * (LY3 * AA + LY4 * Y);

        const mat pnt_pxy = solve(jacob, pnt / 16.);

        for(unsigned J = 0, K = 0, L = 1, M = 2, N = 4; J < s_node; ++J, K += 3, L += 3, M += 3, ++N) {
            m_ip.BP(2, L) = -(m_ip.BP(0, M) = m_ip.BM(0, K) = m_ip.BM(2, L) = pn_pxy(0, J));
            m_ip.BP(1, L) = -(m_ip.BP(2, M) = m_ip.BM(2, K) = m_ip.BM(1, L) = pn_pxy(1, J));
            m_ip.BM(0, M) = pnt_pxy(0, J);
            m_ip.BM(1, M) = pnt_pxy(1, N);
            m_ip.BM(2, M) = pnt_pxy(0, N) + pnt_pxy(1, J);
        }

        auto& s_ip = m_ip.sec_int_pt;
        s_ip.clear();
        s_ip.reserve(t_plan.n_rows);
        for(unsigned J = 0; J < t_plan.n_rows; ++J) {
            const auto t_eccentricity = .5 * t_plan(J, 0) * thickness;
            s_ip.emplace_back(t_eccentricity, .5 * thickness * t_plan(J, 1) * m_plan(I, 2) * det_jacob, mat_proto->get_copy());
            m_stiffness += s_ip.back().factor * m_ip.BM.t() * mat_stiff * m_ip.BM;
            p_stiffness += t_eccentricity * t_eccentricity * s_ip.back().factor * m_ip.BP.t() * mat_stiff * m_ip.BP;
        }
    }

    transform_from_local_to_global(initial_stiffness = reshuffle(m_stiffness, p_stiffness));
    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int S4::update_status() {
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
    mat::fixed<12, 12> m_stiffness(fill::zeros);
    vec3 t_stress;
    vec::fixed<12> m_resistance(fill::zeros);

    auto p_stiffness = penalty_stiffness;
    vec p_resistance = p_stiffness * p_disp;

    for(const auto& I : int_pt) {
        const vec m_strain = I.BM * m_disp, p_strain = I.BP * p_disp;

        for(const auto& J : I.sec_int_pt) if(SUANPAN_SUCCESS != J.s_material->update_trial_status(m_strain + J.eccentricity * p_strain)) return SUANPAN_FAIL;

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

    trial_resistance = reshuffle(m_resistance, p_resistance);
    trial_stiffness = reshuffle(m_stiffness, p_stiffness);

    if(is_nlgeom()) trial_geometry = transform_to_global_geometry(trial_stiffness, trial_resistance, g_disp);

    transform_from_local_to_global(trial_resistance);
    transform_from_local_to_global(trial_stiffness);

    return SUANPAN_SUCCESS;
}

int S4::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.s_material->commit_status();
    return code;
}

int S4::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.s_material->clear_status();
    return code;
}

int S4::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.s_material->reset_status();
    return code;
}

vector<vec> S4::record(const OutputType P) {
    vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) append_to(data, J.s_material->record(P));
    return data;
}

void S4::print() {
    suanpan_info("A four-node planar shell element using GQ12 for membrane action and Mindlin for plate action.\n");
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void S4::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < s_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

void S4::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, s_node, fill::zeros);

    if(OutputType::A == type) t_disp = reshape(get_current_acceleration(), s_dof, s_node);
    else if(OutputType::V == type) t_disp = reshape(get_current_velocity(), s_dof, s_node);
    else if(OutputType::U == type) t_disp = reshape(get_current_displacement(), s_dof, s_node);

    for(unsigned I = 0; I < s_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void S4::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat t_mat = amplifier * mat(reshape(get_current_displacement(), s_dof, s_node)).rows(0, 2).t() + get_coordinate(3);
    for(unsigned I = 0; I < s_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), t_mat(I, 0), t_mat(I, 1), t_mat(I, 2));
}

#endif
