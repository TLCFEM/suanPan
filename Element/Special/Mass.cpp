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

vtkSmartPointer<vtkCell> MassBase::Setup(const uvec& encoding) const {
    auto cell = vtkSmartPointer<vtkVertex>::New();
    for(unsigned I = 0; I < get_node_number(); ++I) cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(encoding(I)));
    return cell;
}

mat MassBase::GetData(const OutputType P) {
    const auto n_dof = get_dof_number();
    const auto n_node = get_node_number();

    if(OutputType::A == P) return resize(reshape(get_current_acceleration(), n_dof, n_node), 6, n_node);
    if(OutputType::V == P) return resize(reshape(get_current_velocity(), n_dof, n_node), 6, n_node);
    if(OutputType::U == P) return resize(reshape(get_current_displacement(), n_dof, n_node), 6, n_node);

    return {};
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
