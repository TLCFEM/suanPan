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

#include "DCP4.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>
#include <Toolbox/utility.h>

const uvec DCP4::u_dof{0, 1, 3, 4, 6, 7, 9, 10};
const uvec DCP4::d_dof{2, 5, 8, 11};

DCP4::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& N, mat&& PNPXY)
    : coor(std::move(C))
    , weight(W)
    , m_material(std::move(M))
    , n_mat(std::move(N))
    , pn_mat(std::move(PNPXY))
    , b_mat(3, 8, fill::zeros) {}

DCP4::DCP4(const unsigned T, uvec&& N, const unsigned M, const double CL, const double RR, const double TH)
    : MaterialElement2D(T, m_node, m_dof, std::move(N), uvec{M}, false, {DOF::U1, DOF::U2, DOF::DMG})
    , release_rate(RR)
    , thickness(TH) { access::rw(characteristic_length) = CL; }

int DCP4::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(PlaneType::E == material_proto->get_plane_type()) suanpan::hacker(thickness) = 1.;

    const auto ele_coor = get_coordinate(2);

    if(0. >= characteristic_length) access::rw(characteristic_length) = 2. * sqrt(area::shoelace(ele_coor));

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

    initial_stiffness.zeros(m_size, m_size);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1)};
        const auto pn = shape::quad(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob), material_proto->get_copy(), shape::quad(t_vec, 0), solve(jacob, pn));

        auto& c_pt = int_pt.back();

        for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += 2, L += 2) {
            c_pt.b_mat(0, K) = c_pt.b_mat(2, L) = c_pt.pn_mat(0, J);
            c_pt.b_mat(2, K) = c_pt.b_mat(1, L) = c_pt.pn_mat(1, J);
        }
        initial_stiffness(u_dof, u_dof) += c_pt.weight * thickness * c_pt.b_mat.t() * ini_stiffness * c_pt.b_mat;
        initial_stiffness(d_dof, d_dof) += c_pt.weight * thickness * release_rate / characteristic_length * c_pt.n_mat.t() * c_pt.n_mat;
        initial_stiffness(d_dof, d_dof) += c_pt.weight * thickness * release_rate * characteristic_length * c_pt.pn_mat.t() * c_pt.pn_mat;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    if(const auto t_density = material_proto->get_density(); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(const auto& I : int_pt) {
            const auto t_factor = t_density * I.weight * thickness;
            for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = J, M = L; K < m_node; ++K, M += m_dof) initial_mass(L, M) += t_factor * I.n_mat(J) * I.n_mat(K);
        }
        for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    return SUANPAN_SUCCESS;
}

int DCP4::update_status() {
    const auto t_disp = get_trial_displacement();

    const vec t_damage = t_disp(d_dof);

    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    for(const auto& I : int_pt) {
        vec t_strain(3, fill::zeros);
        for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
            t_strain(0) += t_disp(K) * I.pn_mat(0, J);
            t_strain(1) += t_disp(L) * I.pn_mat(1, J);
            t_strain(2) += t_disp(K) * I.pn_mat(1, J) + t_disp(L) * I.pn_mat(0, J);
        }

        if(I.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        const auto pow_term = 1. - dot(t_damage, I.n_mat);
        const auto damage = pow(pow_term, 2.);
        const auto t_factor = I.weight * thickness;

        trial_stiffness(u_dof, u_dof) += t_factor * damage * I.b_mat.t() * I.m_material->get_trial_stiffness() * I.b_mat;
        trial_stiffness(u_dof, d_dof) -= t_factor * 2. * pow_term * I.b_mat.t() * I.m_material->get_trial_stress() * I.n_mat;
        trial_stiffness(d_dof, d_dof) += t_factor * (2. * I.maximum_energy + release_rate / characteristic_length) * I.n_mat.t() * I.n_mat;
        trial_stiffness(d_dof, d_dof) += t_factor * release_rate * characteristic_length * I.pn_mat.t() * I.pn_mat;

        trial_resistance(u_dof) += t_factor * damage * I.b_mat.t() * I.m_material->get_trial_stress();
        trial_resistance(d_dof) -= t_factor * 2. * I.n_mat.t() * I.maximum_energy;
    }

    trial_resistance(d_dof) += trial_stiffness(d_dof, d_dof) * t_damage;

    return SUANPAN_SUCCESS;
}

int DCP4::commit_status() {
    auto code = 0;
    for(auto& I : int_pt) {
        I.commit_status(I.m_material);
        I.maximum_energy = std::max(I.maximum_energy, I.strain_energy);
        code += I.m_material->commit_status();
    }
    return code;
}

int DCP4::clear_status() {
    auto code = 0;
    for(auto& I : int_pt) {
        I.clear_status();
        I.maximum_energy = 0.;
        code += I.m_material->clear_status();
    }
    return code;
}

int DCP4::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

vector<vec> DCP4::record(const OutputType P) {
    if(P == OutputType::DAMAGE) return {get_current_displacement()(d_dof)};

    vector<vec> data;
    for(const auto& I : int_pt) append_to(data, I.m_material->record(P));
    return data;
}

void DCP4::print() {
    suanpan_info("A four-node membrane element (DCP4){}.\n", nlgeom ? " with nonlinear geometry (TL formulation)" : "");
    suanpan_info("The nodes connected are:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("IP {}:\t", I + 1);
        suanpan_info(int_pt[I].coor);
        int_pt[I].m_material->print();
    }
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void DCP4::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void DCP4::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration()(u_dof), 2, m_node);
    else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity()(u_dof), 2, m_node);
    else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement()(u_dof), 2, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

mat DCP4::GetData(const OutputType P) {
    if(P == OutputType::DAMAGE) {
        mat t_damage(6, m_node, fill::zeros);
        t_damage.row(0) = get_current_displacement()(d_dof).t();
        return t_damage;
    }

    mat A(int_pt.size(), 4);
    mat B(6, int_pt.size(), fill::zeros);

    for(size_t I = 0; I < int_pt.size(); ++I) {
        if(const auto C = int_pt[I].m_material->record(P); !C.empty()) B(0, I, size(C[0])) = C[0];
        A.row(I) = interpolation::linear(int_pt[I].coor);
    }

    mat data(m_node, 4);

    data.row(0) = interpolation::linear(-1., -1.);
    data.row(1) = interpolation::linear(1., -1.);
    data.row(2) = interpolation::linear(1., 1.);
    data.row(3) = interpolation::linear(-1., 1.);

    return (data * solve(A, B.t())).t();
}

void DCP4::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement()(u_dof), 2, m_node).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
