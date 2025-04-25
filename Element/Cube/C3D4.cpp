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

#include "C3D4.h"
#include <Domain/DomainBase.h>
#include <Material/Material3D/Material3D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/shape.h>
#include <Toolbox/tensor.h>

C3D4::C3D4(const unsigned T, uvec&& N, const unsigned M, const bool F)
    : MaterialElement3D(T, c_node, c_dof, std::move(N), uvec{M}, F) {}

int C3D4::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    c_material = material_proto->get_copy();

    mat ele_coor(c_node, c_node);
    ele_coor.col(0).fill(1.);
    ele_coor.cols(1, 3) = get_coordinate(3);

    access::rw(volume) = det(ele_coor) / 6.;

    const mat inv_coor = inv(ele_coor);

    pn_pxyz = inv_coor.rows(1, 3);

    strain_mat.zeros(6, c_size);
    for(unsigned J = 0, K = 0, L = 1, M = 2; J < c_node; ++J, K += c_dof, L += c_dof, M += c_dof) {
        strain_mat(0, K) = strain_mat(3, L) = strain_mat(5, M) = pn_pxyz(0, J);
        strain_mat(3, K) = strain_mat(1, L) = strain_mat(4, M) = pn_pxyz(1, J);
        strain_mat(5, K) = strain_mat(4, L) = strain_mat(2, M) = pn_pxyz(2, J);
    }

    trial_stiffness = current_stiffness = initial_stiffness = volume * strain_mat.t() * c_material->get_initial_stiffness() * strain_mat;

    const rowvec n = mean(ele_coor) * inv_coor;

    if(const auto t_density = c_material->get_density() * volume; t_density > 0.) {
        initial_mass.zeros(c_size, c_size);
        for(auto I = 0u, K = 0u; I < c_node; ++I, K += c_dof) for(auto J = I, L = K; J < c_node; ++J, L += c_dof) initial_mass(K, L) += t_density * n(I) * n(J);
        for(auto I = 0u, K = 1u; I < c_size; I += c_dof, K += c_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + c_dof, L = K + c_dof; J < c_size; J += c_dof, L += c_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    body_force.zeros(c_size, c_dof);
    for(auto J = 0u, L = 0u; J < c_node; ++J, L += c_dof) for(auto K = 0llu; K < c_dof; ++K) body_force(L + K, K) = n(J);

    return SUANPAN_SUCCESS;
}

int C3D4::update_status() {
    if(nlgeom) {
        const mat ele_disp = reshape(get_trial_displacement(), c_dof, c_node);

        mat BN(6, c_size);

        const mat gradient = ele_disp * pn_pxyz.t() + eye(c_dof, c_dof);
        for(unsigned I = 0, J = 0, K = 1, L = 2; I < c_node; ++I, J += c_dof, K += c_dof, L += c_dof) {
            BN(0, J) = pn_pxyz(0, I) * gradient(0, 0);
            BN(1, J) = pn_pxyz(1, I) * gradient(0, 1);
            BN(2, J) = pn_pxyz(2, I) * gradient(0, 2);
            BN(0, K) = pn_pxyz(0, I) * gradient(1, 0);
            BN(1, K) = pn_pxyz(1, I) * gradient(1, 1);
            BN(2, K) = pn_pxyz(2, I) * gradient(1, 2);
            BN(0, L) = pn_pxyz(0, I) * gradient(2, 0);
            BN(1, L) = pn_pxyz(1, I) * gradient(2, 1);
            BN(2, L) = pn_pxyz(2, I) * gradient(2, 2);
            BN(3, J) = pn_pxyz(0, I) * gradient(0, 1) + pn_pxyz(1, I) * gradient(0, 0);
            BN(4, J) = pn_pxyz(1, I) * gradient(0, 2) + pn_pxyz(2, I) * gradient(0, 1);
            BN(5, J) = pn_pxyz(2, I) * gradient(0, 0) + pn_pxyz(0, I) * gradient(0, 2);
            BN(3, K) = pn_pxyz(0, I) * gradient(1, 1) + pn_pxyz(1, I) * gradient(1, 0);
            BN(4, K) = pn_pxyz(1, I) * gradient(1, 2) + pn_pxyz(2, I) * gradient(1, 1);
            BN(5, K) = pn_pxyz(2, I) * gradient(1, 0) + pn_pxyz(0, I) * gradient(1, 2);
            BN(3, L) = pn_pxyz(0, I) * gradient(2, 1) + pn_pxyz(1, I) * gradient(2, 0);
            BN(4, L) = pn_pxyz(1, I) * gradient(2, 2) + pn_pxyz(2, I) * gradient(2, 1);
            BN(5, L) = pn_pxyz(2, I) * gradient(2, 0) + pn_pxyz(0, I) * gradient(2, 2);
        }

        if(c_material->update_trial_status(tensor::strain::to_voigt(tensor::strain::to_green(gradient))) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        auto& t_stress = c_material->get_trial_stress();

        trial_stiffness = volume * BN.t() * c_material->get_trial_stiffness() * BN;
        trial_resistance = volume * BN.t() * t_stress;

        trial_geometry.zeros(c_size, c_size);
        const auto sigma = tensor::stress::to_tensor(t_stress);

        for(unsigned I = 0, K = 0, M = 1, N = 2; I < c_node; ++I, K += c_dof, M += c_dof, N += c_dof) {
            const vec t_vec = sigma * pn_pxyz.col(I);
            auto t_factor = volume * dot(pn_pxyz.col(I), t_vec);
            trial_geometry(K, K) += t_factor;
            trial_geometry(M, M) += t_factor;
            trial_geometry(N, N) += t_factor;
            for(auto J = I + 1; J < c_node; ++J) {
                t_factor = volume * dot(pn_pxyz.col(J), t_vec);
                auto L = c_dof * J;
                trial_geometry(L, K) += t_factor;
                trial_geometry(K, L) += t_factor;
                trial_geometry(++L, M) += t_factor;
                trial_geometry(M, L) += t_factor;
                trial_geometry(++L, N) += t_factor;
                trial_geometry(N, L) += t_factor;
            }
        }
    }
    else {
        if(c_material->update_trial_status(strain_mat * get_trial_displacement()) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        trial_stiffness = volume * strain_mat.t() * c_material->get_trial_stiffness() * strain_mat;
        trial_resistance = volume * strain_mat.t() * c_material->get_trial_stress();
    }

    return SUANPAN_SUCCESS;
}

int C3D4::commit_status() { return c_material->commit_status(); }

int C3D4::clear_status() { return c_material->clear_status(); }

int C3D4::reset_status() { return c_material->reset_status(); }

std::vector<vec> C3D4::record(const OutputType P) { return c_material->record(P); }

void C3D4::print() {
    suanpan_info("C3D4 element connects:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    c_material->print();
    suanpan_info("Strain:\t", c_material->get_trial_strain());
    suanpan_info("Stress:\t", c_material->get_trial_stress());
}

#ifdef SUANPAN_VTK
#include <vtkTetra.h>

void C3D4::Setup() {
    vtk_cell = vtkSmartPointer<vtkTetra>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < c_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

void C3D4::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, c_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 2) = reshape(get_current_acceleration(), c_dof, c_node);
    else if(OutputType::V == type) t_disp.rows(0, 2) = reshape(get_current_velocity(), c_dof, c_node);
    else if(OutputType::U == type) t_disp.rows(0, 2) = reshape(get_current_displacement(), c_dof, c_node);

    for(unsigned I = 0; I < c_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void C3D4::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(3) + amplifier * reshape(get_current_displacement(), c_dof, c_node).t();
    for(unsigned I = 0; I < c_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
