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

#include "Material3D.h"
#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

Material3D::Material3D(const unsigned T, const double R)
    : Material(T, MaterialType::D3, R) {}

vector<vec> Material3D::record(const OutputType P) {
    if(P == OutputType::S11) return {vec{current_stress(0)}};
    if(P == OutputType::S22) return {vec{current_stress(1)}};
    if(P == OutputType::S33) return {vec{current_stress(2)}};
    if(P == OutputType::S12) return {vec{current_stress(3)}};
    if(P == OutputType::S23) return {vec{current_stress(4)}};
    if(P == OutputType::S13) return {vec{current_stress(5)}};
    if(P == OutputType::E11) return {vec{current_strain(0)}};
    if(P == OutputType::E22) return {vec{current_strain(1)}};
    if(P == OutputType::E33) return {vec{current_strain(2)}};
    if(P == OutputType::E12) return {vec{current_strain(3)}};
    if(P == OutputType::E23) return {vec{current_strain(4)}};
    if(P == OutputType::E13) return {vec{current_strain(5)}};

    if(P == OutputType::SP) {
        vec principal_stress;
        eig_sym(principal_stress, tensor::stress::to_tensor(current_stress));
        return {principal_stress};
    }
    if(P == OutputType::SP1) {
        vec principal_stress;
        eig_sym(principal_stress, tensor::stress::to_tensor(current_stress));
        return {vec{principal_stress(0)}};
    }
    if(P == OutputType::SP2) {
        vec principal_stress;
        eig_sym(principal_stress, tensor::stress::to_tensor(current_stress));
        return {vec{principal_stress(1)}};
    }
    if(P == OutputType::SP3) {
        vec principal_stress;
        eig_sym(principal_stress, tensor::stress::to_tensor(current_stress));
        return {vec{principal_stress(2)}};
    }

    if(P == OutputType::EP) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(current_strain));
        return {principal_strain};
    }
    if(P == OutputType::EP1) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(current_strain));
        return {vec{principal_strain(0)}};
    }
    if(P == OutputType::EP2) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(current_strain));
        return {vec{principal_strain(1)}};
    }
    if(P == OutputType::EP3) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(current_strain));
        return {vec{principal_strain(2)}};
    }

    if(P == OutputType::EEP) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(solve(initial_stiffness, current_stress)));
        return {principal_strain};
    }
    if(P == OutputType::EEP1) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(solve(initial_stiffness, current_stress)));
        return {vec{principal_strain(0)}};
    }
    if(P == OutputType::EEP2) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(solve(initial_stiffness, current_stress)));
        return {vec{principal_strain(1)}};
    }
    if(P == OutputType::EEP3) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(solve(initial_stiffness, current_stress)));
        return {vec{principal_strain(2)}};
    }

    if(P == OutputType::PEP) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(current_strain - solve(initial_stiffness, current_stress)));
        return {principal_strain};
    }
    if(P == OutputType::PEP1) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(current_strain - solve(initial_stiffness, current_stress)));
        return {vec{principal_strain(0)}};
    }
    if(P == OutputType::PEP2) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(current_strain - solve(initial_stiffness, current_stress)));
        return {vec{principal_strain(1)}};
    }
    if(P == OutputType::PEP3) {
        vec principal_strain;
        eig_sym(principal_strain, tensor::strain::to_tensor(current_strain - solve(initial_stiffness, current_stress)));
        return {vec{principal_strain(2)}};
    }

    if(P == OutputType::EEQ) return {vec{sqrt(2. / 3.) * tensor::strain::norm(current_strain)}};
    if(P == OutputType::EEEQ) return {vec{sqrt(2. / 3.) * tensor::strain::norm(solve(initial_stiffness, current_stress))}};
    if(P == OutputType::PEEQ) return {vec{sqrt(2. / 3.) * tensor::strain::norm(current_strain - solve(initial_stiffness, current_stress))}};
    if(P == OutputType::HYDRO) return {vec{tensor::mean3(current_stress)}};
    if(P == OutputType::MISES) return {vec{sqrt(1.5) * tensor::stress::norm(tensor::dev(current_stress))}};

    return Material::record(P);
}
