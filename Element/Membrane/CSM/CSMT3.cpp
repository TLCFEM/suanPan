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

#include "CSMT3.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/utility.h>

const uvec CSMT3::t_dof{0, 1, 3, 4, 6, 7};
const uvec CSMT3::r_dof{2, 5, 8};

CSMT3::IntegrationPoint::IntegrationPoint(rowvec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , weight(W)
    , m_material(std::move(M)) {}

CSMT3::CSMT3(const unsigned T, uvec&& NT, const unsigned MT, const double TH, const double L)
    : MaterialElement2D(T, m_node, m_dof, std::move(NT), uvec{MT}, false, {DOF::U1, DOF::U2, DOF::UR3})
    , thickness(TH) { access::rw(characteristic_length) = L; }

int CSMT3::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(!material_proto->is_support_couple()) {
        suanpan_warning("Element {} is assigned with a material that does not support couple stress.\n", get_tag());
        return SUANPAN_FAIL;
    }

    if(PlaneType::E == material_proto->get_plane_type()) suanpan::hacker(thickness) = 1.;

    mat ele_coor(m_node, m_node);
    ele_coor.col(0).fill(1.);
    ele_coor.cols(1, 2) = get_coordinate(2);

    access::rw(area) = .5 * det(ele_coor);
    if(characteristic_length < 0.) access::rw(characteristic_length) = sqrt(area);

    const mat inv_coor = inv(ele_coor);
    // const mat pn_pxy = inv_coor.rows(1, 2);
    /* partial derivatives of shape functions
     * should be constant for constant strain triangle
     * [ pn1px pn2px pn3px ]
     * [ pn1py pn2py pn3py ]
     */

    const auto& t_size = t_dof.n_elem;
    const auto& r_size = r_dof.n_elem;

    mat l_p(3, t_size, fill::zeros), j_p(1, t_size, fill::zeros), j_q(2, r_size, fill::zeros);

    const auto& j_s = j_p;

    for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += 2, L += 2) {
        // \mathbold{L}\mathbold{\phi}_\mathbold{u}
        // l_p(0, K) = l_p(2, L) = pn_pxy(0, J);
        // l_p(2, K) = l_p(1, L) = pn_pxy(1, J);

        // \mathbold{J}\mathbold{\phi}_\mathbold{u}
        // j_p(0, K) = -pn_pxy(1, J);
        // j_p(0, L) = pn_pxy(0, J);

        // \mathbold{J}\mathbold{\phi}_\mathbold{\theta}
        // j_q(0, J) = pn_pxy(1, J);
        // j_q(1, J) = -pn_pxy(0, J);

        // \mathbold{J}\mathbold{\phi}_\mathbold{\mu}
        // j_s(0, K) = -pn_pxy(1, J);
        // j_s(0, L) = pn_pxy(0, J);

        j_q(1, J) = -(j_p(0, L) = l_p(0, K) = l_p(2, L) = inv_coor(1, J));
        j_p(0, K) = -(j_q(0, J) = l_p(2, K) = l_p(1, L) = inv_coor(2, J));
    }

    j_p *= .5;
    j_q *= .5;

    const mat phi_a = eye(3, 3), phi_b = inv(material_proto->get_initial_stiffness());

    const auto t_factor = area * thickness;

    const mat E2 = t_factor * phi_b.t() * material_proto->get_initial_stiffness() * phi_b;
    const mat H1 = -2. * t_factor * j_p.t() * j_s;
    const mat H2 = t_factor * l_p.t() * phi_a;
    const mat H5 = t_factor * phi_b.t() * phi_a;

    int_pt.clear();
    int_pt.reserve(3);
    int_pt.emplace_back(mean(ele_coor.rows(uvec{1, 2})), t_factor / 3., material_proto->get_copy());
    int_pt.emplace_back(mean(ele_coor.rows(uvec{2, 0})), t_factor / 3., material_proto->get_copy());
    int_pt.emplace_back(mean(ele_coor.rows(uvec{0, 1})), t_factor / 3., material_proto->get_copy());

    mat E1(t_size, t_size, fill::zeros), H3(r_size, t_size, fill::zeros), H4(t_size, t_size, fill::zeros);

    for(auto& I : int_pt) {
        I.m_material->set_characteristic_length(characteristic_length);
        I.m_material->initialize_couple(D);

        const rowvec n = I.coor * inv_coor;

        mat phi_s(2, t_size, fill::zeros);

        const auto& phi_q = n;
        const auto& phi_r = phi_s;

        for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += 2, L += 2)
            // \mathbold{\phi}_\mathbold{\theta}
            // phi_q(0, J) = n(J);

            // \mathbold{\phi}_\mathbold{\mu}
            // phi_s(0, K) = n(J);
            // phi_s(1, L) = n(J);

            phi_s(0, K) = phi_s(1, L) = n(J);

        E1 += I.weight * phi_r.t() * I.m_material->get_initial_couple_stiffness() * phi_r;
        H3 += 2. * I.weight * (phi_q.t() * j_s - j_q.t() * phi_s);
        H4 += I.weight * phi_r.t() * phi_s;

        I.b1 = phi_b;
        I.b2 = I.b3 = phi_r;
    }

    const mat T1 = solve(H5.t(), H2.t());
    const mat T2 = solve(H4.t(), H1.t());
    const mat T3 = solve(H4.t(), H3.t());

    initial_stiffness.set_size(m_size, m_size);
    initial_stiffness(t_dof, t_dof) = T1.t() * E2 * T1 + T2.t() * E1 * T2;
    initial_stiffness(t_dof, r_dof) = T2.t() * E1 * T3;
    initial_stiffness(r_dof, t_dof) = T3.t() * E1 * T2;
    initial_stiffness(r_dof, r_dof) = T3.t() * E1 * T3;

    trial_stiffness = current_stiffness = initial_stiffness;

    for(auto& I : int_pt) {
        I.b1 *= T1;
        I.b2 *= T2;
        I.b3 *= T3;
    }

    return SUANPAN_SUCCESS;
}

int CSMT3::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    for(const auto& I : int_pt) {
        if(SUANPAN_SUCCESS != I.m_material->update_trial_status(I.b1 * t_disp(t_dof))) return SUANPAN_FAIL;
        if(SUANPAN_SUCCESS != I.m_material->update_couple_trial_status(I.b2 * t_disp(t_dof) + I.b3 * t_disp(r_dof))) return SUANPAN_FAIL;

        trial_stiffness(t_dof, t_dof) += I.weight * I.b1.t() * I.m_material->get_trial_stiffness() * I.b1 + I.weight * I.b2.t() * I.m_material->get_trial_couple_stiffness() * I.b2;
        trial_stiffness(t_dof, r_dof) += I.weight * I.b2.t() * I.m_material->get_trial_couple_stiffness() * I.b3;
        trial_stiffness(r_dof, t_dof) += I.weight * I.b3.t() * I.m_material->get_trial_couple_stiffness() * I.b2;
        trial_stiffness(r_dof, r_dof) += I.weight * I.b3.t() * I.m_material->get_trial_couple_stiffness() * I.b3;

        trial_resistance(t_dof) += I.weight * I.b1.t() * I.m_material->get_trial_stress() + I.weight * I.b2.t() * I.m_material->get_trial_couple_stress();
        trial_resistance(r_dof) += I.weight * I.b3.t() * I.m_material->get_trial_couple_stress();
    }

    return SUANPAN_SUCCESS;
}

int CSMT3::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status() + I.m_material->clear_couple_status();
    return code;
}

int CSMT3::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status() + I.m_material->clear_couple_status();
    return code;
}

int CSMT3::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status() + I.m_material->reset_couple_status();
    return code;
}

vector<vec> CSMT3::record(const OutputType P) {
    vector<vec> data;
    for(const auto& I : int_pt) append_to(data, I.m_material->record(P));
    return data;
}

void CSMT3::print() {
    suanpan_info("CSMT3 element connects:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(const auto& I : int_pt) I.m_material->print();
}

#ifdef SUANPAN_VTK
#include <vtkTriangle.h>

void CSMT3::Setup() {
    vtk_cell = vtkSmartPointer<vtkTriangle>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void CSMT3::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration(), m_dof, m_node).eval().head_rows(2);
    else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity(), m_dof, m_node).eval().head_rows(2);
    else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node).eval().head_rows(2);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void CSMT3::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).eval().head_rows(2).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
