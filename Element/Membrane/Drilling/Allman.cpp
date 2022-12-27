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

#include "Allman.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/utility.h>

Allman::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::forward<vec>(C))
    , weight(W)
    , m_material(std::forward<unique_ptr<Material>>(M))
    , strain_mat(3, m_size) {}

mat Allman::form_coor(const mat& C) {
    const auto &X1 = C(0, 0), &X2 = C(1, 0), &X3 = C(2, 0), &Y1 = C(0, 1), &Y2 = C(1, 1), &Y3 = C(2, 1);

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

field<mat> Allman::form_transform(const mat& C) {
    const auto &X1 = C(0, 0), &X2 = C(1, 0), &X3 = C(2, 0), &Y1 = C(0, 1), &Y2 = C(1, 1), &Y3 = C(2, 1);

    mat BMX(6, 9, fill::zeros), BMY(6, 9, fill::zeros);

    BMX(0, 0) = BMX(1, 3) = BMX(2, 6) = 1.;
    BMX(3, 0) = BMX(3, 3) = BMX(4, 3) = BMX(4, 6) = BMX(5, 6) = BMX(5, 0) = .5;
    BMX(3, 2) = -(BMX(3, 5) = .125 * (Y2 - Y1));
    BMX(4, 5) = -(BMX(4, 8) = .125 * (Y3 - Y2));
    BMX(5, 8) = -(BMX(5, 2) = .125 * (Y1 - Y3));

    BMY(0, 1) = BMY(1, 4) = BMY(2, 7) = 1.;
    BMY(3, 1) = BMY(3, 4) = BMY(4, 4) = BMY(4, 7) = BMY(5, 7) = BMY(5, 1) = .5;
    BMY(3, 5) = -(BMY(3, 2) = .125 * (X2 - X1));
    BMY(4, 8) = -(BMY(4, 5) = .125 * (X3 - X2));
    BMY(5, 2) = -(BMY(5, 8) = .125 * (X1 - X3));

    return {BMX, BMY};
}

Allman::Allman(const unsigned T, uvec&& NT, const unsigned MT, const double TH)
    : MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(NT), uvec{MT}, false, {DOF::U1, DOF::U2, DOF::UR3})
    , thickness(TH) {}

int Allman::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get<Material>(material_tag(0));

    if(PlaneType::E == static_cast<PlaneType>(mat_proto->get_parameter(ParameterType::PLANETYPE))) suanpan::hacker(thickness) = 1.;

    auto& mat_stiff = mat_proto->get_initial_stiffness();

    const auto coor = get_coordinate(2);
    const auto ele_coor = form_coor(coor);

    const mat inv_coor = inv(ele_coor);

    area = .5 * det(ele_coor(span(0, 2), span(0, 2)));

    const auto dkt_trans = form_transform(coor);

    const auto &BMX = dkt_trans(0), &BMY = dkt_trans(1);

    initial_stiffness.zeros(m_size, m_size);

    int_pt.clear();
    int_pt.reserve(3);
    for(uword I = 0; I < 3; ++I) {
        int_pt.emplace_back(vec{ele_coor(I + 3, 1), ele_coor(I + 3, 2)}, area * thickness / 3., mat_proto->get_copy());

        auto& c_pt = int_pt.back();

        const mat pn_pxy = shape::triangle(c_pt.coor, 1) * inv_coor;
        c_pt.strain_mat.row(0) = pn_pxy.row(0) * BMX;
        c_pt.strain_mat.row(1) = pn_pxy.row(1) * BMY;
        c_pt.strain_mat.row(2) = pn_pxy.row(0) * BMY + pn_pxy.row(1) * BMX;

        initial_stiffness += c_pt.weight * c_pt.strain_mat.t() * mat_stiff * c_pt.strain_mat;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    if(const auto t_density = mat_proto->get_parameter(ParameterType::DENSITY); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(const auto& I : int_pt) {
            const rowvec n_int = shape::triangle(I.coor, 0) * inv_coor;
            const auto t_factor = t_density * I.weight;
            for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = J, P = L; K < m_node; ++K, P += m_dof) initial_mass(L, P) += t_factor * n_int(J) * n_int(K);
        }
        for(unsigned I = 0, K = 1; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    return SUANPAN_SUCCESS;
}

int Allman::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    for(const auto& I : int_pt) {
        if(I.m_material->update_trial_status(I.strain_mat * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        trial_stiffness += I.weight * I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat;
        trial_resistance += I.weight * I.strain_mat.t() * I.m_material->get_trial_stress();
    }

    return SUANPAN_SUCCESS;
}

int Allman::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int Allman::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int Allman::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

vector<vec> Allman::record(const OutputType T) {
    vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.m_material->record(T)) data.emplace_back(J);
    return data;
}

void Allman::print() { node_encoding.t().print("Allman element connects:"); }

#ifdef SUANPAN_VTK
#include <vtkTriangle.h>

void Allman::Setup() {
    vtk_cell = vtkSmartPointer<vtkTriangle>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void Allman::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_acceleration(), m_dof, m_node);
    else if(OutputType::V == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_velocity(), m_dof, m_node);
    else if(OutputType::U == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_displacement(), m_dof, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void Allman::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * mat(reshape(get_current_displacement(), m_dof, m_node)).rows(0, 1).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
