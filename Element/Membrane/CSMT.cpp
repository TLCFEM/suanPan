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

#include "CSMT.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/utility.h>
#include <Toolbox/tensorToolbox.h>

CSMT::CSMT(const unsigned T, uvec&& NT, const unsigned MT, const double TH)
	: MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(NT), uvec{MT}, false)
	, thickness(TH) {}

void CSMT::initialize(const shared_ptr<DomainBase>& D) {
	auto& material_proto = D->get<Material>(material_tag(0));

	if(suanpan::approx_equal(static_cast<double>(PlaneType::E), material_proto->get_parameter(ParameterType::PLANETYPE))) suanpan::hacker(thickness) = 1.;

	m_material = material_proto->get_copy();

	mat ele_coor(m_node, m_node);
	ele_coor.col(0).fill(1.);
	ele_coor.cols(1, 2) = get_coordinate(2);

	access::rw(area) = .5 * det(ele_coor);
	access::rw(characteristic_length) = sqrt(area);

	m_material->set_characteristic_length(characteristic_length);
	m_material->initialize_couple(D);

	const mat inv_coor = inv(ele_coor);
	const rowvec n = mean(ele_coor) * inv_coor;
	const mat pn_pxy = inv_coor.rows(1, 2);
	/* partial derivatives of shape functions
	 * [ pn1px pn2px pn3px ]
	 * [ pn1py pn2py pn3py ]
	 */

	mat phi_q(1, 3, fill::zeros), phi_s(2, 6, fill::zeros);
	mat l_p(3, 6, fill::zeros), j_p(1, 6, fill::zeros), j_q(2, 3, fill::zeros);

	const mat phi_r = eye(2, 2), phi_a = eye(3, 3), phi_b = inv(m_material->get_initial_stiffness());

	for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += 2, L += 2) {
		// \mathbold{L}\mathbold{\phi}_\mathbold{u}
		l_p(0, K) = l_p(2, L) = pn_pxy(0, J);
		l_p(2, K) = l_p(1, L) = pn_pxy(1, J);

		// \mathbold{J}\mathbold{\phi}_\mathbold{u}
		j_p(0, K) = -pn_pxy(1, J);
		j_p(0, L) = pn_pxy(0, J);

		// \mathbold{\phi}_\mathbold{\theta}
		phi_q(0, J) = n(J);

		// \mathbold{J}\mathbold{\phi}_\mathbold{\theta}
		j_q(0, J) = pn_pxy(1, J);
		j_q(1, J) = -pn_pxy(0, J);

		// \mathbold{\phi}_\mathbold{\mu}
		phi_s(0, K) = n(J);
		phi_s(1, L) = n(J);

		// \mathbold{J}\mathbold{\phi}_\mathbold{\mu}
		// j_s(0, K) = -pn_pxy(1, J);
		// j_s(0, L) = pn_pxy(0, J);
	}

	j_p *= .5;
	j_q *= .5;

	const auto& j_s = j_p;

	const auto t_factor = area * thickness;

	const mat E1 = t_factor * phi_r.t() * m_material->get_initial_couple_stiffness() * phi_r;
	const mat E2 = t_factor * phi_b.t() * m_material->get_initial_stiffness() * phi_b;
	const mat H1 = -2. * t_factor * j_p.t() * j_s;
	const mat H2 = t_factor * l_p.t() * phi_a;
	const mat H3 = 2. * t_factor * (phi_q.t() * j_s - j_q.t() * phi_s);
	const mat H4 = -2. * t_factor * phi_r.t() * phi_s;
	const mat H5 = t_factor * phi_b.t() * phi_a;

	const mat T1 = solve(H5.t(), H2.t());
	const mat T2 = solve(H4.t(), H1.t());
	const mat T3 = solve(H4.t(), H3.t());

	initial_stiffness.set_size(m_size, m_size);
	initial_stiffness(t_dof, t_dof) = T1.t() * E2 * T1 + T2.t() * E1 * T2;
	initial_stiffness(t_dof, r_dof) = T2.t() * E1 * T3;
	initial_stiffness(r_dof, t_dof) = T3.t() * E1 * T2;
	initial_stiffness(r_dof, r_dof) = T3.t() * E1 * T3;

	trial_stiffness = current_stiffness = initial_stiffness;

	if(const auto t_density = area * thickness * m_material->get_parameter(ParameterType::DENSITY); t_density > 0.) {
		initial_mass.zeros(m_size, m_size);
		for(auto I = 0u, K = 0u; I < m_node; ++I, K += m_dof) for(auto J = I, L = K; J < m_node; ++J, L += m_dof) initial_mass(K, L) += t_density * n(I) * n(J);
		for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
			initial_mass(K, K) = initial_mass(I, I);
			for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
		}
		ConstantMass(this);
	}
}

int CSMT::update_status() {
	trial_resistance = trial_stiffness * get_trial_displacement();

	return SUANPAN_SUCCESS;
}

int CSMT::commit_status() { return m_material->commit_status(); }

int CSMT::clear_status() { return m_material->clear_status(); }

int CSMT::reset_status() { return m_material->reset_status(); }

vector<vec> CSMT::record(const OutputType T) { return m_material->record(T); }

void CSMT::print() {
	suanpan_info("CSMT element.\n");
	if(!is_initialized()) return;
	suanpan_info("Material model response:\n");
	m_material->print();
}

#ifdef SUANPAN_VTK
#include <vtkTriangle.h>

void CSMT::Setup() {
	vtk_cell = vtkSmartPointer<vtkTriangle>::New();
	const auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < m_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void CSMT::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
	mat t_disp(6, m_node, fill::zeros);

	if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration(), m_dof, m_node);
	else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity(), m_dof, m_node);
	else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);

	for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

mat CSMT::GetData(const OutputType P) {
	vec t_stress(6, fill::zeros);
	if(const auto t_data = m_material->record(P); !t_data.empty()) t_stress(uvec{0, 1, 3}) = t_data[0];
	return repmat(t_stress, 1, m_node);
}

void CSMT::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).eval().head_rows(2).t();
	for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
