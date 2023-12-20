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

#include "PCPE4DC.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>
#include <Toolbox/tensor.h>
#include <Toolbox/utility.h>

const uvec PCPE4DC::s_dof{0, 1, 4, 5, 8, 9, 12, 13};
const uvec PCPE4DC::f_dof{2, 3, 6, 7, 10, 11, 14, 15};

PCPE4DC::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , weight(W)
    , m_material(std::move(M))
    , strain_mat(3, 2llu * m_node, fill::zeros) {}

PCPE4DC::PCPE4DC(const unsigned T, uvec&& N, const unsigned MS, const unsigned MF, const double AL, const double NN, const double KK)
    : MaterialElement2D(T, m_node, m_dof, std::move(N), uvec{MS, MF}, false)
    , alpha(AL)
    , porosity(NN)
    , k(KK) {}

int PCPE4DC::initialize(const shared_ptr<DomainBase>& D) {
    auto& s_mat = D->get<Material>(material_tag(0));
    auto& f_mat = D->get<Material>(material_tag(1));

    // validate material type
    if(PlaneType::E != s_mat->get_plane_type()) {
        suanpan_error("Only plane strain material for solid phase is supported.\n");
        return SUANPAN_FAIL;
    }
    if(MaterialType::DS != f_mat->get_material_type()) {
        suanpan_error("Only isotropic fluid phase is supported.\n");
        return SUANPAN_FAIL;
    }

    // validate bulk modulus
    const auto ks = s_mat->get_parameter(ParameterType::BULKMODULUS);
    const auto kf = f_mat->get_parameter(ParameterType::BULKMODULUS);

    if(suanpan::approx_equal(ks, 0.) || suanpan::approx_equal(kf, 0.)) {
        suanpan_error("A zero bulk modulus is detected.\n");
        return SUANPAN_FAIL;
    }

    q = ks * kf / (porosity * ks + (alpha - porosity) * kf);

    // compute density ratio
    const auto s_density = (1. - porosity) * s_mat->get_density();
    const auto f_density = porosity * f_mat->get_density();
    const auto s_ratio = s_density / (s_density + f_density);
    const auto f_ratio = f_density / (s_density + f_density);

    const auto ele_coor = get_coordinate(2);

    access::rw(characteristic_length) = sqrt(area::shoelace(ele_coor));

    auto& ini_stiffness = s_mat->get_initial_stiffness();

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

    initial_stiffness.zeros(m_size, m_size);
    body_force.zeros(m_size, 2);

    mat meta_a(m_node, m_node, fill::zeros);
    mat meta_b(2llu * m_node, 2llu * m_node, fill::zeros);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(auto I = 0u; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1)};
        const auto n = compute_shape_function(t_vec, 0);
        const auto pn = compute_shape_function(t_vec, 1);
        const mat jacob = pn * ele_coor;
        const mat pn_pxy = solve(jacob, pn);
        int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob), s_mat->get_copy());

        auto& c_pt = int_pt.back();

        for(auto J = 0u; J < m_node; ++J) {
            const auto K = 2 * J, L = K + 1, M = 4 * J, N = M + 1, O = N + 1, P = O + 1;
            c_pt.strain_mat(0, K) = c_pt.strain_mat(2, L) = pn_pxy(0, J);
            c_pt.strain_mat(2, K) = c_pt.strain_mat(1, L) = pn_pxy(1, J);
            body_force(M, 0) += c_pt.weight * s_ratio * n(J);
            body_force(N, 1) += c_pt.weight * s_ratio * n(J);
            body_force(O, 0) += c_pt.weight * f_ratio * n(J);
            body_force(P, 1) += c_pt.weight * f_ratio * n(J);
        }

        initial_stiffness(s_dof, s_dof) += c_pt.weight * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;
        meta_a += c_pt.weight * n.t() * n;
        meta_b += c_pt.weight * c_pt.strain_mat.head_rows(2).t() * ones(2, 2) * c_pt.strain_mat.head_rows(2);
    }
    meta_k.set_size(m_size, m_size);
    meta_k(s_dof, s_dof) = q * pow(alpha - porosity, 2.) * meta_b;
    meta_k(s_dof, f_dof) = q * porosity * (alpha - porosity) * meta_b;
    meta_k(f_dof, s_dof) = q * porosity * (alpha - porosity) * meta_b;
    meta_k(f_dof, f_dof) = q * porosity * porosity * meta_b;
    trial_stiffness = current_stiffness = initial_stiffness += meta_k;

    const uvec s_dof_a{0, 4, 8, 12};
    const uvec s_dof_b{1, 5, 9, 13};
    const uvec f_dof_a{2, 6, 10, 14};
    const uvec f_dof_b{3, 7, 11, 15};

    initial_viscous.zeros(m_size, m_size);
    initial_viscous(s_dof_a, s_dof_a) = porosity * porosity / k * meta_a;
    initial_viscous(s_dof_b, s_dof_b) = initial_viscous(s_dof_a, s_dof_a);
    initial_viscous(f_dof_a, f_dof_a) = initial_viscous(s_dof_a, s_dof_a);
    initial_viscous(f_dof_b, f_dof_b) = initial_viscous(s_dof_a, s_dof_a);
    initial_viscous(s_dof_a, f_dof_a) = -initial_viscous(s_dof_a, s_dof_a);
    initial_viscous(s_dof_b, f_dof_b) = -initial_viscous(s_dof_a, s_dof_a);
    initial_viscous(f_dof_a, s_dof_a) = -initial_viscous(s_dof_a, s_dof_a);
    initial_viscous(f_dof_b, s_dof_b) = -initial_viscous(s_dof_a, s_dof_a);
    ConstantDamping(this);

    initial_mass.zeros(m_size, m_size);
    initial_mass(s_dof_a, s_dof_a) = s_density * meta_a;
    initial_mass(s_dof_b, s_dof_b) = initial_mass(s_dof_a, s_dof_a);
    initial_mass(f_dof_a, f_dof_a) = f_density * meta_a;
    initial_mass(f_dof_b, f_dof_b) = initial_mass(f_dof_a, f_dof_a);
    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

int PCPE4DC::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_stiffness = meta_k;
    trial_resistance = meta_k * t_disp;

    for(const auto& I : int_pt) {
        if(I.m_material->update_trial_status(I.strain_mat * t_disp(s_dof)) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_stiffness(s_dof, s_dof) += I.weight * I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat;
        trial_resistance(s_dof) += I.weight * I.strain_mat.t() * I.m_material->get_trial_stress();
    }

    trial_viscous_force = trial_viscous * get_trial_velocity();

    return SUANPAN_SUCCESS;
}

int PCPE4DC::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int PCPE4DC::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int PCPE4DC::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

mat PCPE4DC::compute_shape_function(const mat& coordinate, const unsigned order) const { return shape::quad(coordinate, order, m_node); }

vector<vec> PCPE4DC::record(const OutputType P) {
    vector<vec> data;

    if(P == OutputType::PP) {
        const auto t_disp = get_current_displacement();
        for(const auto& I : int_pt) data.emplace_back(vec{q * tensor::trace2(I.strain_mat * ((porosity - alpha) * t_disp(s_dof) - porosity * t_disp(f_dof)))});
    }
    else for(const auto& I : int_pt) append_to(data, I.m_material->record(P));

    return data;
}

void PCPE4DC::print() {
    suanpan_info("A four-node membrane element (PCPE4DC).\n");
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

void PCPE4DC::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void PCPE4DC::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration()(s_dof), 2, m_node);
    else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity()(s_dof), 2, m_node);
    else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement()(s_dof), 2, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

mat PCPE4DC::GetData(const OutputType P) {
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

void PCPE4DC::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement()(s_dof), 2, m_node).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
