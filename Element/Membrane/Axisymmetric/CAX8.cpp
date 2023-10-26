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

#include "CAX8.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>

CAX8::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::forward<vec>(C))
    , weight(W)
    , m_material(std::forward<unique_ptr<Material>>(M))
    , strain_mat(4, m_size, fill::zeros) {}

vec CAX8::isoparametric_mapping(const vec& in) {
    const auto &X1 = in(0), &X2 = in(1), &X3 = in(2), &X4 = in(3), &X5 = in(4), &X6 = in(5), &X7 = in(6), &X8 = in(7);

    vec out(8);

    out(0) = .25 * (2. * (X5 + X6 + X7 + X8) - X1 - X2 - X3 - X4);
    out(1) = .5 * (X6 - X8);
    out(2) = .5 * (X7 - X5);
    out(3) = .25 * (X1 + X2 + X3 + X4 - 2. * (X5 + X7));
    out(4) = .25 * (X1 - X2 + X3 - X4);
    out(5) = .25 * (X1 + X2 + X3 + X4 - 2. * (X6 + X8));
    out(6) = .25 * (2. * (X5 - X7) - X1 - X2 + X3 + X4);
    out(7) = .25 * (2. * (X8 - X6) - X1 + X2 + X3 - X4);

    return out;
}

CAX8::CAX8(const unsigned T, uvec&& N, const unsigned M, const bool R, const bool F)
    : MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, F)
    , reduced_scheme(R) {}

int CAX8::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(PlaneType::A != material_proto->get_plane_type()) {
        suanpan_warning("Element {} is assigned with an inconsistent material.\n", get_tag());
        return SUANPAN_FAIL;
    }

    const auto ele_coor = get_coordinate(m_dof);

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const IntegrationPlan plan(2, reduced_scheme ? 2 : 3, IntegrationType::GAUSS);

    initial_stiffness.zeros(m_size, m_size);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        const auto& X = plan(I, 0);
        const auto& Y = plan(I, 1);
        vec t_vec{X, Y};
        const auto n = shape::quad(t_vec, 0, m_node);
        const auto pn = shape::quad(t_vec, 1, m_node);
        const mat jacob = pn * ele_coor;
        const mat pn_pxy = solve(jacob, pn);
        const auto gx = dot(vec{1., X, Y, X * X, X * Y, Y * Y, X * X * Y, X * Y * Y}, isoparametric_mapping(ele_coor.col(0)));

        int_pt.emplace_back(std::move(t_vec), 2. * datum::pi * gx * plan(I, 2) * det(jacob), material_proto->get_copy());

        auto& c_int_pt = int_pt.back();

        for(auto J = 0u, K = 0u, L = 1u; J < m_node; ++J, K += m_dof, L += m_dof) {
            c_int_pt.strain_mat(0, K) = c_int_pt.strain_mat(3, L) = pn_pxy(0, J);
            c_int_pt.strain_mat(3, K) = c_int_pt.strain_mat(1, L) = pn_pxy(1, J);
            c_int_pt.strain_mat(2, K) = n(J) / gx;
        }
        initial_stiffness += c_int_pt.weight * c_int_pt.strain_mat.t() * ini_stiffness * c_int_pt.strain_mat;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    // trial_mass = current_mass = initial_mass.zeros(m_size, m_size);

    return SUANPAN_SUCCESS;
}

int CAX8::update_status() {
    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    for(const auto& I : int_pt) {
        vec t_strain(4, fill::zeros);
        for(unsigned J = 0, K = 0; J < m_node; ++J, K += m_dof) {
            const auto& t_disp = node_ptr[J].lock()->get_trial_displacement();
            t_strain(0) += t_disp(0) * I.strain_mat(0, K);
            t_strain(1) += t_disp(1) * I.strain_mat(3, K);
            t_strain(2) += t_disp(0) * I.strain_mat(2, K);
            t_strain(3) += t_disp(0) * I.strain_mat(3, K) + t_disp(1) * I.strain_mat(0, K);
        }
        if(I.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_stiffness += I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat * I.weight;
        trial_resistance += I.strain_mat.t() * I.m_material->get_trial_stress() * I.weight;
    }

    return SUANPAN_SUCCESS;
}

int CAX8::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int CAX8::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int CAX8::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

vector<vec> CAX8::record(const OutputType T) {
    vector<vec> data;
    for(const auto& I : int_pt) append_to(data, I.m_material->record(T));
    return data;
}

void CAX8::print() {
    suanpan_info("A CAX8{} element{}.\n", reduced_scheme ? "R" : "", nlgeom ? " with nonlinear geometry on" : "");
    suanpan_info("The nodes connected are:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(const auto& I : int_pt) I.m_material->print();
}

#ifdef SUANPAN_VTK
#include <vtkQuadraticQuad.h>

void CAX8::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuadraticQuad>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void CAX8::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration(), m_dof, m_node);
    else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity(), m_dof, m_node);
    else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void CAX8::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
