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

#include "C3D20.h"
#include <Domain/DomainBase.h>
#include <Material/Material3D/Material3D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shape.h>
#include <Toolbox/tensor.h>

C3D20::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& P)
    : coor(std::move(C))
    , weight(W)
    , c_material(std::move(M))
    , pn_pxyz(std::move(P))
    , strain_mat(6, c_size) {
    for(auto I = 0u, J = 0u, K = 1u, L = 2u; I < c_node; ++I, J += c_dof, K += c_dof, L += c_dof) {
        strain_mat(0, J) = strain_mat(3, K) = strain_mat(5, L) = pn_pxyz(0, I);
        strain_mat(3, J) = strain_mat(1, K) = strain_mat(4, L) = pn_pxyz(1, I);
        strain_mat(5, J) = strain_mat(4, K) = strain_mat(2, L) = pn_pxyz(2, I);
    }
}

C3D20::C3D20(const unsigned T, uvec&& N, const unsigned M, const bool R, const bool F)
    : MaterialElement3D(T, c_node, c_dof, std::move(N), uvec{M}, F)
    , reduced_scheme(R) {}

int C3D20::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    const auto ele_coor = get_coordinate(c_dof);

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const IntegrationPlan plan(3, reduced_scheme ? 2 : 3, reduced_scheme ? IntegrationType::GAUSS : IntegrationType::IRONS);

    initial_stiffness.zeros(c_size, c_size);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1), plan(I, 2)};
        const auto pn = compute_shape_function(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(std::move(t_vec), plan(I, 3) * det(jacob), material_proto->get_copy(), solve(jacob, pn));

        const auto& c_pt = int_pt.back();
        initial_stiffness += c_pt.weight * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    if(const auto t_density = material_proto->get_density(); t_density > 0.) {
        initial_mass.zeros(c_size, c_size);
        for(const auto& I : int_pt) {
            const auto n_int = compute_shape_function(I.coor, 0);
            const auto t_factor = t_density * I.weight;
            for(auto J = 0u, L = 0u; J < c_node; ++J, L += c_dof) for(auto K = J, M = L; K < c_node; ++K, M += c_dof) initial_mass(L, M) += t_factor * n_int(J) * n_int(K);
        }
        for(auto I = 0u, K = 1u, L = 2u; I < c_size; I += c_dof, K += c_dof, L += c_dof) {
            initial_mass(K, K) = initial_mass(L, L) = initial_mass(I, I);
            for(auto J = I + c_dof, M = J + 1, N = J + 2; J < c_size; J += c_dof, M += c_dof, N += c_dof) initial_mass(J, I) = initial_mass(K, M) = initial_mass(L, N) = initial_mass(M, K) = initial_mass(N, L) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    body_force.zeros(c_size, c_dof);
    for(const auto& I : int_pt) {
        const mat n_int = I.weight * compute_shape_function(I.coor, 0);
        for(auto J = 0u, L = 0u; J < c_node; ++J, L += c_dof) for(auto K = 0llu; K < c_dof; ++K) body_force(L + K, K) += n_int(J);
    }

    return SUANPAN_SUCCESS;
}

int C3D20::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_stiffness.zeros(c_size, c_size);
    trial_resistance.zeros(c_size);

    if(nlgeom) {
        const mat ele_disp = reshape(t_disp, c_dof, c_node);

        trial_geometry.zeros(c_size, c_size);

        mat BN(6, c_size);
        for(const auto& I : int_pt) {
            const mat gradient = ele_disp * I.pn_pxyz.t() + eye(c_dof, c_dof);
            for(unsigned J = 0, K = 0, L = 1, M = 2; J < c_node; ++J, K += c_dof, L += c_dof, M += c_dof) {
                BN(0, K) = I.pn_pxyz(0, J) * gradient(0, 0);
                BN(1, K) = I.pn_pxyz(1, J) * gradient(0, 1);
                BN(2, K) = I.pn_pxyz(2, J) * gradient(0, 2);
                BN(0, L) = I.pn_pxyz(0, J) * gradient(1, 0);
                BN(1, L) = I.pn_pxyz(1, J) * gradient(1, 1);
                BN(2, L) = I.pn_pxyz(2, J) * gradient(1, 2);
                BN(0, M) = I.pn_pxyz(0, J) * gradient(2, 0);
                BN(1, M) = I.pn_pxyz(1, J) * gradient(2, 1);
                BN(2, M) = I.pn_pxyz(2, J) * gradient(2, 2);
                BN(3, K) = I.pn_pxyz(0, J) * gradient(0, 1) + I.pn_pxyz(1, J) * gradient(0, 0);
                BN(4, K) = I.pn_pxyz(1, J) * gradient(0, 2) + I.pn_pxyz(2, J) * gradient(0, 1);
                BN(5, K) = I.pn_pxyz(2, J) * gradient(0, 0) + I.pn_pxyz(0, J) * gradient(0, 2);
                BN(3, L) = I.pn_pxyz(0, J) * gradient(1, 1) + I.pn_pxyz(1, J) * gradient(1, 0);
                BN(4, L) = I.pn_pxyz(1, J) * gradient(1, 2) + I.pn_pxyz(2, J) * gradient(1, 1);
                BN(5, L) = I.pn_pxyz(2, J) * gradient(1, 0) + I.pn_pxyz(0, J) * gradient(1, 2);
                BN(3, M) = I.pn_pxyz(0, J) * gradient(2, 1) + I.pn_pxyz(1, J) * gradient(2, 0);
                BN(4, M) = I.pn_pxyz(1, J) * gradient(2, 2) + I.pn_pxyz(2, J) * gradient(2, 1);
                BN(5, M) = I.pn_pxyz(2, J) * gradient(2, 0) + I.pn_pxyz(0, J) * gradient(2, 2);
            }

            if(I.c_material->update_trial_status(tensor::strain::to_voigt(tensor::strain::to_green(gradient))) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

            auto& t_stress = I.c_material->get_trial_stress();

            const auto sigma = tensor::stress::to_tensor(t_stress);

            for(unsigned J = 0, K = 0, L = 1, M = 2; J < c_node; ++J, K += c_dof, L += c_dof, M += c_dof) {
                const vec t_vec = I.weight * sigma * I.pn_pxyz.col(J);
                auto t_factor = dot(I.pn_pxyz.col(J), t_vec);
                trial_geometry(K, K) += t_factor;
                trial_geometry(L, L) += t_factor;
                trial_geometry(M, M) += t_factor;
                for(auto N = J + 1, O = c_dof * N, P = O + 1, Q = P + 1; N < c_node; ++N, O += c_dof, P += c_dof, Q += c_dof) {
                    t_factor = dot(I.pn_pxyz.col(N), t_vec);
                    trial_geometry(O, K) += t_factor;
                    trial_geometry(P, L) += t_factor;
                    trial_geometry(Q, M) += t_factor;
                    trial_geometry(K, O) += t_factor;
                    trial_geometry(L, P) += t_factor;
                    trial_geometry(M, Q) += t_factor;
                }
            }

            trial_stiffness += I.weight * BN.t() * I.c_material->get_trial_stiffness() * BN;
            trial_resistance += I.weight * BN.t() * t_stress;
        }
    }
    else
        for(const auto& I : int_pt) {
            if(I.c_material->update_trial_status(I.strain_mat * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
            trial_stiffness += I.weight * I.strain_mat.t() * I.c_material->get_trial_stiffness() * I.strain_mat;
            trial_resistance += I.weight * I.strain_mat.t() * I.c_material->get_trial_stress();
        }

    return SUANPAN_SUCCESS;
}

int C3D20::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->commit_status();
    return code;
}

int C3D20::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->clear_status();
    return code;
}

int C3D20::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->reset_status();
    return code;
}

mat C3D20::compute_shape_function(const mat& coordinate, const unsigned order) const { return shape::cube(coordinate, order, c_node); }

vector<vec> C3D20::record(const OutputType P) {
    vector<vec> data;
    for(const auto& I : int_pt) append_to(data, I.c_material->record(P));
    return data;
}

void C3D20::print() {
    suanpan_info("A C3D20 element{}{}.\n", reduced_scheme ? " reduced eight-point Gauss integration" : " full fourteen-point Irons integration", nlgeom ? " nonlinear geometry" : "");
    suanpan_info("The element connects nodes:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(const auto& t_pt : int_pt) {
        t_pt.c_material->print();
        suanpan_info("Strain:\t", t_pt.c_material->get_trial_strain());
        suanpan_info("Stress:\t", t_pt.c_material->get_trial_stress());
    }
}

#ifdef SUANPAN_VTK
#include <vtkQuadraticHexahedron.h>

void C3D20::Setup() {
    vtk_cell = vtkSmartPointer<vtkQuadraticHexahedron>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < c_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

void C3D20::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, c_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 2) = reshape(get_current_acceleration(), c_dof, c_node);
    else if(OutputType::V == type) t_disp.rows(0, 2) = reshape(get_current_velocity(), c_dof, c_node);
    else if(OutputType::U == type) t_disp.rows(0, 2) = reshape(get_current_displacement(), c_dof, c_node);

    for(unsigned I = 0; I < c_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void C3D20::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(3) + amplifier * reshape(get_current_displacement(), c_dof, c_node).t();
    for(unsigned I = 0; I < c_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
