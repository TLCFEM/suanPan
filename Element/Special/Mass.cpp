/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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

Mass::Mass(const unsigned T, const unsigned NT, const double MA, uvec&& DT)
	: Element(T, 1, static_cast<unsigned>(DT.max()), uvec{NT})
	, magnitude(MA)
	, dof_label(DT - 1) {}

void Mass::initialize(const shared_ptr<DomainBase>&) {
	initial_mass.zeros(dof_label.max() + 1, dof_label.max() + 1);
	for(const auto& I : dof_label) initial_mass(I, I) = magnitude;

	ConstantMass(this);
}

int Mass::update_status() { return SUANPAN_SUCCESS; }

int Mass::commit_status() { return SUANPAN_SUCCESS; }

int Mass::clear_status() { return SUANPAN_SUCCESS; }

int Mass::reset_status() { return SUANPAN_SUCCESS; }

void Mass::print() { suanpan_info("A point mass element.\n"); }

#ifdef SUANPAN_VTK
#include <vtkVertex.h>

void Mass::Setup() {
	vtk_cell = vtkSmartPointer<vtkVertex>::New();
	const auto ele_coor = get_coordinate(3);
	for(unsigned I = 0; I < get_node_number(); ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(0, 2));
	}
}

void Mass::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
	const auto n_dof = get_dof_number();
	const auto n_node = get_node_number();

	mat t_disp(6, n_node, fill::zeros);

	if(OutputType::A == type) t_disp.rows(0, n_dof - 1) = reshape(get_current_acceleration(), n_dof, n_node);
	else if(OutputType::V == type) t_disp.rows(0, n_dof - 1) = reshape(get_current_velocity(), n_dof, n_node);
	else if(OutputType::U == type) t_disp.rows(0, n_dof - 1) = reshape(get_current_displacement(), n_dof, n_node);

	for(unsigned I = 0; I < n_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

mat Mass::GetData(const OutputType) { return {}; }

void Mass::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const auto n_dof = get_dof_number();
	const auto n_node = get_node_number();

	mat ele_disp = get_coordinate(n_dof) + amplifier * reshape(get_current_displacement(), n_dof, n_node).t();
	ele_disp.resize(n_node, 3);
	for(unsigned I = 0; I < n_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
