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

#include "PCPE4UC.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>

PCPE4UC::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::forward<vec>(C))
    , weight(W)
    , m_material(std::forward<unique_ptr<Material>>(M))
    , strain_mat(3, m_size, fill::zeros) {}

PCPE4UC::PCPE4UC(const unsigned T, uvec&& N, const unsigned MS, const unsigned MF, const double AL, const double NN)
    : MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(N), uvec{MS, MF}, false)
    , alpha(AL)
    , porosity(NN) {}

int PCPE4UC::initialize(const shared_ptr<DomainBase>& D) {
    auto& s_mat = D->get<Material>(material_tag(0));
    auto& f_mat = D->get<Material>(material_tag(1));

    // validate material type
    if(PlaneType::E != static_cast<PlaneType>(s_mat->get_parameter(ParameterType::PLANETYPE))) {
        suanpan_error("PCPE4UC %u only supports the plane strain material for solid phase.\n", get_tag());
        return SUANPAN_FAIL;
    }
    if(MaterialType::DS != f_mat->get_material_type()) {
        suanpan_error("PCPE4UC %u only supports the isotropic fluid phase.\n", get_tag());
        return SUANPAN_FAIL;
    }

    // validate bulk modulus
    const auto ks = s_mat->get_parameter(ParameterType::BULKMODULUS);
    const auto kf = f_mat->get_parameter(ParameterType::BULKMODULUS);

    if(suanpan::approx_equal(ks, 0.)) {
        suanpan_error("PCPE4UC %u solid phase returns a zero bulk modulus.\n", get_tag());
        return SUANPAN_FAIL;
    }
    if(suanpan::approx_equal(kf, 0.)) {
        suanpan_error("PCPE4UC %u fluid phase returns a zero bulk modulus.\n", get_tag());
        return SUANPAN_FAIL;
    }

    q = ks * kf / (porosity * ks + (alpha - porosity) * kf);

    const auto ele_coor = get_coordinate(2);

    access::rw(characteristic_length) = sqrt(area::shoelace(ele_coor));

    auto& ini_stiffness = s_mat->get_initial_stiffness();

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

    meta_k.zeros(m_size, m_size);
    initial_stiffness.zeros(m_size, m_size);
    body_force.zeros(m_size, 2);

    mat meta_a(m_node, m_node, fill::zeros);

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
            const auto K = m_dof * J, L = K + 1;
            c_pt.strain_mat(0, K) = c_pt.strain_mat(2, L) = pn_pxy(0, J);
            c_pt.strain_mat(2, K) = c_pt.strain_mat(1, L) = pn_pxy(1, J);
            body_force(K, 0) += c_pt.weight * n(J);
            body_force(L, 1) += c_pt.weight * n(J);
        }

        initial_stiffness += c_pt.weight * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;
        meta_a += c_pt.weight * n.t() * n;
        meta_k += c_pt.weight * q * alpha * alpha * c_pt.strain_mat.head_rows(2).t() * ones(2, 2) * c_pt.strain_mat.head_rows(2);
    }
    trial_stiffness = current_stiffness = initial_stiffness += meta_k;

    const uvec s_dof_a{0, 2, 4, 6};
    const uvec s_dof_b{1, 3, 5, 7};
    initial_mass.zeros(m_size, m_size);
    initial_mass(s_dof_a, s_dof_a) = ((1. - porosity) * s_mat->get_parameter(ParameterType::DENSITY) + porosity * f_mat->get_parameter(ParameterType::DENSITY)) * meta_a;
    initial_mass(s_dof_b, s_dof_b) = initial_mass(s_dof_a, s_dof_a);
    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

int PCPE4UC::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_stiffness = meta_k;
    trial_resistance = meta_k * t_disp;

    for(const auto& I : int_pt) {
        if(I.m_material->update_trial_status(I.strain_mat * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_stiffness += I.weight * I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat;
        trial_resistance += I.weight * I.strain_mat.t() * I.m_material->get_trial_stress();
    }

    return SUANPAN_SUCCESS;
}

int PCPE4UC::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int PCPE4UC::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int PCPE4UC::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

mat PCPE4UC::compute_shape_function(const mat& coordinate, const unsigned order) const { return shape::quad(coordinate, order, m_node); }

vector<vec> PCPE4UC::record(const OutputType P) {
    vector<vec> output;

    if(P == OutputType::NMISES) {
        mat A(int_pt.size(), 4);
        vec B(int_pt.size(), fill::zeros);

        for(size_t I = 0; I < int_pt.size(); ++I) {
            if(const auto C = int_pt[I].m_material->record(OutputType::MISES); !C.empty()) B(I) = C.cbegin()->at(0);
            A.row(I) = interpolation::linear(int_pt[I].coor);
        }

        const vec X = solve(A, B);

        output.emplace_back(vec{dot(interpolation::linear(-1., -1.), X)});
        output.emplace_back(vec{dot(interpolation::linear(1., -1.), X)});
        output.emplace_back(vec{dot(interpolation::linear(1., 1.), X)});
        output.emplace_back(vec{dot(interpolation::linear(-1., 1.), X)});
    }
    else if(P == OutputType::PP) {
        const auto t_disp = get_current_displacement();
        for(const auto& I : int_pt) output.emplace_back(vec{-alpha * q * tensor::trace2(I.strain_mat * t_disp)});
    }
    else for(const auto& I : int_pt) for(const auto& J : I.m_material->record(P)) output.emplace_back(J);

    return output;
}

void PCPE4UC::print() {
    suanpan_info("Element %u is a four-node membrane element (PCPE4UC).\n", get_tag());
    node_encoding.t().print("The nodes connected are:");
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("Integration Point %llu:\t", I + 1);
        int_pt[I].coor.t().print();
        int_pt[I].m_material->print();
    }
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void PCPE4UC::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuad>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void PCPE4UC::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration(), 2, m_node);
    else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity(), 2, m_node);
    else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement(), 2, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

mat PCPE4UC::GetData(const OutputType P) {
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

void PCPE4UC::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), 2, m_node).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
