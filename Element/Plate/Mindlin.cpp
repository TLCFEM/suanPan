/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "Mindlin.h"
#include <Domain/DomainBase.h>
#include <Material/Material.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/utility.h>

Mindlin::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const double E, const double F, unique_ptr<Material>&& M)
    : eccentricity(E)
    , factor(F)
    , p_material(std::forward<unique_ptr<Material>>(M)) {}

Mindlin::IntegrationPoint::SectionIntegrationPoint::SectionIntegrationPoint(const SectionIntegrationPoint& old_obj)
    : eccentricity(old_obj.eccentricity)
    , factor(old_obj.factor)
    , p_material(old_obj.p_material->get_copy()) {}

Mindlin::IntegrationPoint::IntegrationPoint(vec&& C)
    : coor(std::forward<vec>(C))
    , strain_mat(3, p_size, fill::zeros) {}

Mindlin::Mindlin(const unsigned T, uvec&& NT, const unsigned MT, const double TH, const unsigned IPN)
    : MaterialElement2D(T, p_node, p_dof, std::forward<uvec>(NT), uvec{MT}, false, {DOF::U1, DOF::U2, DOF::UR3})
    , thickness(TH)
    , num_section_ip(IPN) {}

int Mindlin::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get<Material>(material_tag(0));

    const auto ele_coor = get_coordinate(2);

    auto& mat_stiff = mat_proto->get_initial_stiffness();

    auto shear_modulus = mat_proto->get_parameter(ParameterType::G);
    if(suanpan::approx_equal(0., shear_modulus)) shear_modulus = mat_proto->get_parameter(ParameterType::SHEARMODULUS);
    if(suanpan::approx_equal(0., shear_modulus)) shear_modulus = .5 * mat_proto->get_parameter(ParameterType::E) / (1. + mat_proto->get_parameter(ParameterType::POISSONSRATIO));
    if(suanpan::approx_equal(0., shear_modulus)) shear_modulus = mat_stiff.at(2, 2);
    if(suanpan::approx_equal(0., shear_modulus)) shear_modulus = mat_proto->get_parameter(ParameterType::E);

    // reduced integration for the Kirchhoff constraint
    vec t_vec(2, fill::zeros);
    const auto n = shape::quad(t_vec, 0);
    auto pn = shape::quad(t_vec, 1);
    mat jacob = pn * ele_coor;
    mat pn_pxy = solve(jacob, pn);
    mat penalty_mat(2, p_size, fill::zeros);
    for(uword I = 0; I < p_node; ++I) {
        penalty_mat(0, I * p_dof) = pn_pxy(1, I);
        penalty_mat(1, I * p_dof) = pn_pxy(0, I);
        penalty_mat(0, I * p_dof + 1llu) = -(penalty_mat(1, I * p_dof + 2llu) = n(I));
    }
    initial_stiffness = penalty_stiffness = 10. / 3. * shear_modulus * thickness * det(jacob) * penalty_mat.t() * penalty_mat;

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);
    const IntegrationPlan sec_plan(1, num_section_ip, IntegrationType::GAUSS);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        int_pt.emplace_back(vec{plan(I, 0), plan(I, 1)});

        pn = shape::quad(int_pt.back().coor, 1);
        jacob = pn * ele_coor;
        pn_pxy = solve(jacob, pn);
        const auto det_jacob = det(jacob);

        auto& strain_mat = int_pt.back().strain_mat;
        for(unsigned J = 0, K = 1, L = 2; J < p_node; ++J, K += p_dof, L += p_dof) {
            strain_mat(2, K) = -(strain_mat(0, L) = pn_pxy(0, J));
            strain_mat(1, K) = -(strain_mat(2, L) = pn_pxy(1, J));
        }

        auto& current_ip = int_pt.back().sec_int_pt;
        current_ip.clear();
        current_ip.reserve(num_section_ip);
        for(unsigned J = 0; J < num_section_ip; ++J) {
            const auto t_eccentricity = .5 * sec_plan(J, 0) * thickness;
            current_ip.emplace_back(t_eccentricity, .5 * thickness * sec_plan(J, 1) * plan(I, 2) * det_jacob, mat_proto->get_copy());
            initial_stiffness += t_eccentricity * t_eccentricity * current_ip.back().factor * strain_mat.t() * mat_stiff * strain_mat;
        }
    }

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int Mindlin::update_status() {
    const auto trial_disp = get_trial_displacement();

    trial_resistance = (trial_stiffness = penalty_stiffness) * trial_disp;
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

int Mindlin::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.p_material->commit_status();
    return code;
}

int Mindlin::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.p_material->clear_status();
    return code;
}

int Mindlin::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) code += J.p_material->reset_status();
    return code;
}

vector<vec> Mindlin::record(const OutputType P) {
    vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.sec_int_pt) for(const auto& K : J.p_material->record(P)) data.emplace_back(K);
    return data;
}

void Mindlin::print() { node_encoding.t().print("A Mindlin plate element connects:"); }

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void Mindlin::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < p_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void Mindlin::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, p_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(2, 4) = reshape(get_current_acceleration(), p_dof, p_node);
    else if(OutputType::V == type) t_disp.rows(2, 4) = reshape(get_current_velocity(), p_dof, p_node);
    else if(OutputType::U == type) t_disp.rows(2, 4) = reshape(get_current_displacement(), p_dof, p_node);

    for(unsigned I = 0; I < p_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void Mindlin::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const auto ele_coor = get_coordinate(2);
    const mat ele_disp = reshape(get_current_displacement(), p_dof, p_node);
    for(unsigned I = 0; I < p_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_coor(I, 0), ele_coor(I, 1), amplifier * ele_disp(0, I));
}

#endif
