/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#include "NonlocalCP4.h"

#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>
#include <Toolbox/utility.h>

const uvec NonlocalCP4::u_dof{0, 1, 3, 4, 6, 7, 9, 10};
const uvec NonlocalCP4::d_dof{2, 5, 8, 11};

NonlocalCP4::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& N, const mat& PNPXY)
    : coor(std::move(C))
    , weight(W)
    , m_material(std::move(M))
    , n_mat(std::move(N))
    , b_mat(4, m_size, fill::zeros) {
    for(auto I{0u}, J{0u}, K{1u}, L{2u}; I < m_node; ++I, J += m_dof, K += m_dof, L += m_dof) {
        b_mat(0, J) = b_mat(2, K) = PNPXY(0, I);
        b_mat(2, J) = b_mat(1, K) = PNPXY(1, I);
        b_mat(3, L) = n_mat(I);
    }
}

NonlocalCP4::NonlocalCP4(const unsigned T, uvec&& N, const unsigned M, const double TH)
    : MaterialElement2D(T, m_node, m_dof, std::move(N), uvec{M}, false, {Node::DOF::U1, Node::DOF::U2, Node::DOF::DAMAGE})
    , thickness(TH) {}

int NonlocalCP4::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(material_proto->nonlocal_size() != 1u) {
        suanpan_error("Element {} requires a nonlocal material with size 1, but got {}.\n", get_tag(), material_proto->nonlocal_size());
        return SUANPAN_FAIL;
    }

    if(PlaneType::E == material_proto->get_plane_type()) suanpan::hacker(thickness) = 1.;

    const auto ele_coor = get_coordinate(2);

    access::rw(characteristic_length) = 2. * std::sqrt(area::shoelace(ele_coor));

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const IntegrationPlan plan(2, 2, IntegrationPlan::Type::GAUSS);

    initial_stiffness.zeros(m_size, m_size);

    const_mat.zeros(d_dof.n_elem, d_dof.n_elem);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I{0}; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1)};
        const auto pn = shape::quad(t_vec, 1);
        const mat jacob = pn * ele_coor;
        const mat pn_mat = solve(jacob, pn);
        auto& c_pt = int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob) * thickness, material_proto->unique_copy(), shape::quad(t_vec, 0), pn_mat);

        const_mat -= c_pt.weight * c_pt.n_mat.t() * c_pt.n_mat + c_pt.weight * characteristic_length * characteristic_length * pn_mat.t() * pn_mat;

        initial_stiffness += c_pt.weight * c_pt.b_mat.t() * ini_stiffness * c_pt.b_mat;
    }
    initial_stiffness(d_dof, d_dof) = const_mat;
    trial_stiffness = current_stiffness = initial_stiffness;

    if(const auto t_density = material_proto->get_density(); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(const auto& I : int_pt) {
            const auto t_factor = t_density * I.weight;
            for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof)
                for(auto K = J, M = L; K < m_node; ++K, M += m_dof) initial_mass(L, M) += t_factor * I.n_mat(J) * I.n_mat(K);
        }
        for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    return SUANPAN_SUCCESS;
}

int NonlocalCP4::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    for(const auto& I : int_pt) {
        if(I.m_material->update_trial_status(I.b_mat * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_resistance += I.weight * I.b_mat.t() * I.m_material->get_trial_stress();
        trial_stiffness += I.weight * I.b_mat.t() * I.m_material->get_trial_stiffness() * I.b_mat;
    }

    trial_stiffness(d_dof, d_dof) = const_mat;
    trial_resistance(d_dof) += const_mat * t_disp(d_dof);

    return SUANPAN_SUCCESS;
}

int NonlocalCP4::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int NonlocalCP4::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int NonlocalCP4::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

std::vector<vec> NonlocalCP4::record(const OutputType P) const {
    if(OutputType::DAMAGE == P) return {get_current_displacement()(d_dof)};

    std::vector<vec> data;
    for(const auto& I : int_pt) suanpan::append_to(data, I.m_material->record(P));
    return data;
}

void NonlocalCP4::print() {
    suanpan_info("A four-node membrane element (NonlocalCP4){}.\n", nlgeom ? " with nonlinear geometry (TL formulation)" : "");
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

vtkSmartPointer<vtkCell> NonlocalCP4::GetCell() const { return vtkSmartPointer<vtkQuad>::New(); }

mat NonlocalCP4::GetData(const OutputType P) {
    if(OutputType::A == P) return resize(reshape(get_current_acceleration()(u_dof), 2, m_node), 6, m_node);
    if(OutputType::V == P) return resize(reshape(get_current_velocity()(u_dof), 2, m_node), 6, m_node);
    if(OutputType::U == P) return resize(reshape(get_current_displacement()(u_dof), 2, m_node), 6, m_node);

    if(OutputType::DAMAGE == P) return get_current_displacement()(d_dof).t();

    mat A(static_cast<uword>(int_pt.size()), 4);
    mat B(6, static_cast<uword>(int_pt.size()), fill::zeros);

    for(uword I{0}; I < int_pt.size(); ++I) {
        if(auto C = int_pt[I].m_material->record(P); !C.empty()) B.col(I) = C[0].resize(6);
        A.row(I) = interpolation::linear(int_pt[I].coor);
    }

    mat data(m_node, 4);

    data.row(0) = interpolation::linear(-1., -1.);
    data.row(1) = interpolation::linear(1., -1.);
    data.row(2) = interpolation::linear(1., 1.);
    data.row(3) = interpolation::linear(-1., 1.);

    return (data * solve(A, B.t())).t();
}

mat NonlocalCP4::GetDeformation(const double amplifier) { return get_coordinate(2).t() + amplifier * reshape(get_current_displacement()(u_dof), 2, m_node); }

#endif
