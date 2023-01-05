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

#include "GQ12.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>
#include <Toolbox/utility.h>

GQ12::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::forward<vec>(C))
    , weight(W)
    , m_material(std::forward<unique_ptr<Material>>(M))
    , strain_mat(3, m_size, fill::zeros) {}

GQ12::GQ12(const unsigned T, uvec&& N, const unsigned M, const double TH)
    : MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, false, {DOF::U1, DOF::U2, DOF::UR3})
    , thickness(TH) {}

int GQ12::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get<Material>(material_tag(0));

    if(PlaneType::E == static_cast<PlaneType>(mat_proto->get_parameter(ParameterType::PLANETYPE))) suanpan::hacker(thickness) = 1.;

    auto& mat_stiff = mat_proto->get_initial_stiffness();

    const auto ele_coor = get_coordinate(2);

    access::rw(characteristic_length) = sqrt(area::shoelace(ele_coor));

    const auto LX1 = ele_coor(1, 1) - ele_coor(0, 1);
    const auto LX2 = ele_coor(2, 1) - ele_coor(1, 1);
    const auto LX3 = ele_coor(3, 1) - ele_coor(2, 1);
    const auto LX4 = ele_coor(0, 1) - ele_coor(3, 1);
    const auto LY1 = ele_coor(0, 0) - ele_coor(1, 0);
    const auto LY2 = ele_coor(1, 0) - ele_coor(2, 0);
    const auto LY3 = ele_coor(2, 0) - ele_coor(3, 0);
    const auto LY4 = ele_coor(3, 0) - ele_coor(0, 0);

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

    mat pnt(2, 8);

    initial_stiffness.zeros(m_size, m_size);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1)};
        const auto pn = compute_shape_function(t_vec, 1);
        const mat jacob = pn * ele_coor;
        const mat pn_pxy = solve(jacob, pn);
        int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob) * thickness, mat_proto->get_copy());

        const auto TX = 2. * plan(I, 0);
        const auto TY = 2. * plan(I, 1);

        const auto AA = plan(I, 0) + 1.;
        const auto BB = plan(I, 0) - 1.;
        const auto CC = plan(I, 1) + 1.;
        const auto DD = plan(I, 1) - 1.;

        pnt(0, 0) = +DD * (LX4 * CC - LX1 * TX);
        pnt(0, 1) = +DD * (LX2 * CC + LX1 * TX);
        pnt(0, 2) = -CC * (LX2 * DD - LX3 * TX);
        pnt(0, 3) = -CC * (LX4 * DD + LX3 * TX);
        pnt(0, 4) = +DD * (LY4 * CC - LY1 * TX);
        pnt(0, 5) = +DD * (LY2 * CC + LY1 * TX);
        pnt(0, 6) = -CC * (LY2 * DD - LY3 * TX);
        pnt(0, 7) = -CC * (LY4 * DD + LY3 * TX);
        pnt(1, 0) = -BB * (LX1 * AA - LX4 * TY);
        pnt(1, 1) = +AA * (LX1 * BB + LX2 * TY);
        pnt(1, 2) = +AA * (LX3 * BB - LX2 * TY);
        pnt(1, 3) = -BB * (LX3 * AA + LX4 * TY);
        pnt(1, 4) = -BB * (LY1 * AA - LY4 * TY);
        pnt(1, 5) = +AA * (LY1 * BB + LY2 * TY);
        pnt(1, 6) = +AA * (LY3 * BB - LY2 * TY);
        pnt(1, 7) = -BB * (LY3 * AA + LY4 * TY);

        const mat pnt_pxy = solve(jacob, 625E-4 * pnt);

        auto& c_pt = int_pt.back();

        for(unsigned J = 0, K = 0, L = 1, M = 2, N = 4; J < m_node; ++J, K += m_dof, L += m_dof, M += m_dof, ++N) {
            c_pt.strain_mat(0, K) = c_pt.strain_mat(2, L) = pn_pxy(0, J);
            c_pt.strain_mat(2, K) = c_pt.strain_mat(1, L) = pn_pxy(1, J);
            c_pt.strain_mat(0, M) = pnt_pxy(0, J);
            c_pt.strain_mat(1, M) = pnt_pxy(1, N);
            c_pt.strain_mat(2, M) = pnt_pxy(0, N) + pnt_pxy(1, J);
        }

        initial_stiffness += c_pt.strain_mat.t() * mat_stiff * c_pt.strain_mat * c_pt.weight;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    if(const auto t_density = mat_proto->get_parameter(ParameterType::DENSITY); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(const auto& I : int_pt) {
            const auto n_int = compute_shape_function(I.coor, 0);
            const auto t_factor = t_density * I.weight;
            for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = J, M = L; K < m_node; ++K, M += m_dof) initial_mass(L, M) += t_factor * n_int(J) * n_int(K);
        }
        for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    body_force.zeros(m_size, m_dof);
    for(const auto& I : int_pt) {
        const mat n_int = I.weight * compute_shape_function(I.coor, 0);
        for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = 0llu; K < m_dof; ++K) body_force(L + K, K) += n_int(J);
    }

    return SUANPAN_SUCCESS;
}

int GQ12::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);
    for(const auto& I : int_pt) {
        if(SUANPAN_SUCCESS != I.m_material->update_trial_status(I.strain_mat * t_disp)) return SUANPAN_FAIL;
        trial_stiffness += I.weight * I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat;
        trial_resistance += I.weight * I.strain_mat.t() * I.m_material->get_trial_stress();
    }

    return SUANPAN_SUCCESS;
}

int GQ12::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int GQ12::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int GQ12::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

mat GQ12::compute_shape_function(const mat& coordinate, const unsigned order) const { return shape::quad(coordinate, order, m_node); }

vector<vec> GQ12::record(const OutputType T) {
    vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.m_material->record(T)) data.emplace_back(J);
    return data;
}

void GQ12::print() {
    suanpan_info("A GQ12 element.\n");
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("IP {}:\t", I + 1);
        int_pt[I].m_material->print();
    }
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void GQ12::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void GQ12::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_acceleration(), m_dof, m_node);
    else if(OutputType::V == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_velocity(), m_dof, m_node);
    else if(OutputType::U == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_displacement(), m_dof, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

mat GQ12::GetData(const OutputType P) {
    mat A(int_pt.size(), 4);
    mat B(int_pt.size(), 6, fill::zeros);

    for(size_t I = 0; I < int_pt.size(); ++I) {
        if(const auto C = int_pt[I].m_material->record(P); !C.empty()) B(I, 0, size(C[0])) = C[0];
        A.row(I) = interpolation::linear(int_pt[I].coor);
    }

    mat data(m_node, 4);

    data.row(0) = interpolation::linear(-1., -1.);
    data.row(1) = interpolation::linear(1., -1.);
    data.row(2) = interpolation::linear(1., 1.);
    data.row(3) = interpolation::linear(-1., 1.);

    return (data * solve(A, B)).t();
}

void GQ12::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
