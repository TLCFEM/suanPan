/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "CP6.h"
#include <Domain/DOF.h>
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>

CP6::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& P)
    : coor(std::forward<vec>(C))
    , weight(W)
    , m_material(std::forward<unique_ptr<Material>>(M))
    , pn_pxy(std::forward<mat>(P))
    , strain_mat(3, m_size, fill::zeros) {}

CP6::CP6(const unsigned T, uvec&& NT, const unsigned MT, const double TH, const bool R)
    : MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(NT), uvec{MT}, R, {DOF::U1, DOF::U2})
    , thickness(TH) {}

int CP6::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(PlaneType::E == static_cast<PlaneType>(material_proto->get_parameter(ParameterType::PLANETYPE))) suanpan::hacker(thickness) = 1.;

    mat ele_coor(m_node, m_node, fill::none);
    ele_coor.col(0).fill(1.);
    ele_coor.cols(1, 2) = get_coordinate(2);
    ele_coor.col(3) = ele_coor.col(1) % ele_coor.col(2);
    ele_coor.col(4) = square(ele_coor.col(1));
    ele_coor.col(5) = square(ele_coor.col(2));

    const mat inv_coor = inv(ele_coor);

    area = .5 * det(ele_coor(span(0, 2), span(0, 2)));

    access::rw(characteristic_length) = sqrt(area);

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    initial_stiffness.zeros(m_size, m_size);

    int_pt.clear();
    int_pt.reserve(3);
    for(auto I = 0; I < 3; ++I) {
        vec coor{ele_coor(I + 3llu, 1), ele_coor(I + 3llu, 2)};
        int_pt.emplace_back(std::move(coor), area * thickness / 3., material_proto->get_copy(), shape::triangle(coor, 1) * inv_coor);

        auto& c_pt = int_pt.back();

        for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
            c_pt.strain_mat(0, K) = c_pt.strain_mat(2, L) = c_pt.pn_pxy(0, J);
            c_pt.strain_mat(2, K) = c_pt.strain_mat(1, L) = c_pt.pn_pxy(1, J);
        }
        initial_stiffness += c_pt.weight * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    if(const auto t_density = material_proto->get_parameter(ParameterType::DENSITY); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(const auto& I : int_pt) {
            const rowvec n_int = shape::triangle(I.coor, 0) * inv_coor;
            const auto t_factor = t_density * I.weight;
            for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = J, M = L; K < m_node; ++K, M += m_dof) initial_mass(L, M) += t_factor * n_int(J) * n_int(K);
        }
        for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    body_force.zeros(m_size, m_dof);
    for(const auto& I : int_pt) {
        const rowvec n_int = I.weight * shape::triangle(I.coor, 0) * inv_coor;
        for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = 0llu; K < m_dof; ++K) body_force(L + K, K) += n_int(J);
    }

    return SUANPAN_SUCCESS;
}

int CP6::update_status() {
    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    if(nlgeom) {
        trial_geometry.zeros(m_size, m_size);

        const mat ele_disp = reshape(get_trial_displacement(), m_dof, m_node);

        mat BN(3, m_size);
        for(const auto& I : int_pt) {
            const mat gradient = ele_disp * I.pn_pxy.t() + eye(m_dof, m_dof);
            for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
                BN(0, K) = I.pn_pxy(0, J) * gradient(0, 0);
                BN(1, K) = I.pn_pxy(1, J) * gradient(0, 1);
                BN(0, L) = I.pn_pxy(0, J) * gradient(1, 0);
                BN(1, L) = I.pn_pxy(1, J) * gradient(1, 1);
                BN(2, K) = I.pn_pxy(0, J) * gradient(0, 1) + I.pn_pxy(1, J) * gradient(0, 0);
                BN(2, L) = I.pn_pxy(0, J) * gradient(1, 1) + I.pn_pxy(1, J) * gradient(1, 0);
            }

            if(I.m_material->update_trial_status(tensor::strain::to_voigt(tensor::strain::to_green(gradient))) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

            auto& t_stress = I.m_material->get_trial_stress();

            const auto sigma = tensor::stress::to_tensor(t_stress);

            for(unsigned J = 0, L = 0, M = 1; J < m_node; ++J, L += m_dof, M += m_dof) {
                const vec t_vec = sigma * I.pn_pxy.col(J);
                auto t_factor = I.weight * dot(I.pn_pxy.col(J), t_vec);
                trial_geometry(L, L) += t_factor;
                trial_geometry(M, M) += t_factor;
                for(auto K = J + 1, N = m_dof * K, O = N + 1; K < m_node; ++K, N += m_dof, O += m_dof) {
                    t_factor = I.weight * dot(I.pn_pxy.col(K), t_vec);
                    trial_geometry(N, L) += t_factor;
                    trial_geometry(O, M) += t_factor;
                    trial_geometry(L, N) += t_factor;
                    trial_geometry(M, O) += t_factor;
                }
            }

            trial_stiffness += I.weight * BN.t() * I.m_material->get_trial_stiffness() * BN;
            trial_resistance += I.weight * BN.t() * t_stress;
        }
    }
    else
        for(const auto& I : int_pt) {
            vec t_strain(3, fill::zeros);
            for(unsigned J = 0; J < m_node; ++J) {
                const auto& t_disp = node_ptr[J].lock()->get_trial_displacement();
                t_strain(0) += t_disp(0) * I.pn_pxy(0, J);
                t_strain(1) += t_disp(1) * I.pn_pxy(1, J);
                t_strain(2) += t_disp(0) * I.pn_pxy(1, J) + t_disp(1) * I.pn_pxy(0, J);
            }
            if(I.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

            trial_stiffness += I.weight * I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat;
            trial_resistance += I.weight * I.strain_mat.t() * I.m_material->get_trial_stress();
        }

    return SUANPAN_SUCCESS;
}

int CP6::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int CP6::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int CP6::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

vector<vec> CP6::record(const OutputType T) {
    vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.m_material->record(T)) data.emplace_back(J);
    return data;
}

void CP6::print() {
    node_encoding.t().print("CP6 element connects:");
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(const auto& I : int_pt) I.m_material->print();
}

#ifdef SUANPAN_VTK
#include <vtkQuadraticTriangle.h>

void CP6::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuadraticTriangle>::New();
    const auto ele_coor = get_coordinate(2);
    for(unsigned I = 0; I < m_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
    }
}

void CP6::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, m_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration(), m_dof, m_node);
    else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity(), m_dof, m_node);
    else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);

    for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void CP6::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
    for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
