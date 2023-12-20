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

#include "Mass.h"
#include <Domain/DOF.h>

MassBase::MassBase(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, std::vector<DOF>&& DI)
    : Element(T, NN, ND, std::move(NT), std::move(DI)) { modify_mass = false; }

int MassBase::update_status() { return SUANPAN_SUCCESS; }

int MassBase::commit_status() { return SUANPAN_SUCCESS; }

int MassBase::clear_status() { return SUANPAN_SUCCESS; }

int MassBase::reset_status() { return SUANPAN_SUCCESS; }

void MassBase::print() {
    suanpan_info("A point mass element.\n");
}

#ifdef SUANPAN_VTK
#include <vtkVertex.h>

void MassBase::Setup() {
    vtk_cell = vtkSmartPointer<vtkVertex>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < get_node_number(); ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(0, 2));
    }
}

void MassBase::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    const auto n_dof = get_dof_number();
    const auto n_node = get_node_number();

    mat t_disp(6, n_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, n_dof - 1) = reshape(get_current_acceleration(), n_dof, n_node);
    else if(OutputType::V == type) t_disp.rows(0, n_dof - 1) = reshape(get_current_velocity(), n_dof, n_node);
    else if(OutputType::U == type) t_disp.rows(0, n_dof - 1) = reshape(get_current_displacement(), n_dof, n_node);

    for(unsigned I = 0; I < n_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

mat MassBase::GetData(const OutputType) { return {}; }

void MassBase::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const auto n_dof = get_dof_number();
    const auto n_node = get_node_number();

    mat current_disp(3, n_node, fill::none);
    for(unsigned I = 0; I < n_node; ++I) nodes->GetPoint(static_cast<vtkIdType>(node_encoding(I)), current_disp.colptr(I));

    current_disp.head_rows(n_dof) = get_coordinate(n_dof).t() + amplifier * reshape(get_current_displacement(), n_dof, n_node);
    for(unsigned I = 0; I < n_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), current_disp(0, I), current_disp(1, I), current_disp(2, I));
}

#endif

Mass2D::Mass2D(const unsigned T, const unsigned NT, const double MA, uvec&& DT)
    : MassBase(T, 1, std::min(3u, static_cast<unsigned>(DT.max())), uvec{NT}, [&] {
        std::vector DI(std::min(3llu, DT.max()), DOF::NONE);

        for(const auto I : DT)
            if(1 == I) DI[0] = DOF::U1;
            else if(2 == I) DI[1] = DOF::U2;
            else if(3 == I) DI[2] = DOF::UR3;

        return DI;
    }())
    , magnitude(MA)
    , dof_label(DT - 1) {}

int Mass2D::initialize(const shared_ptr<DomainBase>&) {
    initial_mass.zeros(dof_label.max() + 1, dof_label.max() + 1);
    for(const auto& I : dof_label) initial_mass(I, I) = magnitude;

    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

Mass3D::Mass3D(const unsigned T, const unsigned NT, const double MA, uvec&& DT)
    : MassBase(T, 1, std::min(6u, static_cast<unsigned>(DT.max())), uvec{NT}, [&] {
        std::vector DI(std::min(6llu, DT.max()), DOF::NONE);

        for(const auto I : DT)
            if(1 == I) DI[0] = DOF::U1;
            else if(2 == I) DI[1] = DOF::U2;
            else if(3 == I) DI[2] = DOF::U3;
            else if(4 == I) DI[3] = DOF::UR1;
            else if(5 == I) DI[4] = DOF::UR2;
            else if(6 == I) DI[5] = DOF::UR3;

        return DI;
    }())
    , magnitude(MA)
    , dof_label(DT - 1) {}

int Mass3D::initialize(const shared_ptr<DomainBase>&) {
    initial_mass.zeros(dof_label.max() + 1, dof_label.max() + 1);
    for(const auto& I : dof_label) initial_mass(I, I) = magnitude;

    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

MassPoint2D::MassPoint2D(const unsigned T, const unsigned NT, const double TM)
    : MassBase(T, 1, 2, uvec{NT}, {DOF::U1, DOF::U2})
    , translational_magnitude(TM)
    , rotational_magnitude(0.) {}

MassPoint2D::MassPoint2D(const unsigned T, const unsigned NT, const double TM, const double RM)
    : MassBase(T, 1, 3, uvec{NT}, {DOF::U1, DOF::U2, DOF::UR3})
    , translational_magnitude(TM)
    , rotational_magnitude(RM) {}

int MassPoint2D::initialize(const shared_ptr<DomainBase>&) {
    const auto dof_size = get_dof_number();
    initial_mass.zeros(dof_size, dof_size);
    initial_mass.diag().fill(translational_magnitude);
    if(3u == dof_size) initial_mass(2, 2) = rotational_magnitude;

    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

MassPoint3D::MassPoint3D(const unsigned T, const unsigned NT, const double TM)
    : MassBase(T, 1, 3, uvec{NT}, {DOF::U1, DOF::U2, DOF::U3})
    , translational_magnitude(TM)
    , rotational_magnitude(0.) {}

MassPoint3D::MassPoint3D(const unsigned T, const unsigned NT, const double TM, const double RM)
    : MassBase(T, 1, 6, uvec{NT}, {DOF::U1, DOF::U2, DOF::U3, DOF::UR1, DOF::UR2, DOF::UR3})
    , translational_magnitude(TM)
    , rotational_magnitude(RM) {}

int MassPoint3D::initialize(const shared_ptr<DomainBase>&) {
    const auto dof_size = get_dof_number();
    initial_mass.zeros(dof_size, dof_size);
    initial_mass.diag().fill(translational_magnitude);
    if(6u == dof_size) initial_mass(3, 3) = initial_mass(4, 4) = initial_mass(5, 5) = rotational_magnitude;

    ConstantMass(this);

    return SUANPAN_SUCCESS;
}
