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

#include "Material3D.h"
#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

Material3D::Material3D(const unsigned T, const double R)
    : Material(T, MaterialType::D3, R) {}

vector<vec> Material3D::record(const OutputType P) {
    vector<vec> data;

    if(P == OutputType::S) data.emplace_back(current_stress);
    else if(P == OutputType::HYDRO) data.emplace_back(vec{tensor::mean3(current_stress)});
    else if(P == OutputType::MISES) data.emplace_back(vec{sqrt(1.5) * tensor::stress::norm(tensor::dev(current_stress))});
    else if(P == OutputType::E) data.emplace_back(current_strain);
    else if(P == OutputType::EE) data.emplace_back(solve(initial_stiffness, current_stress));
    else if(P == OutputType::PE) data.emplace_back(current_strain - solve(initial_stiffness, current_stress));
    else if(P == OutputType::EEQ) data.emplace_back(vec{sqrt(2. / 3.) * tensor::strain::norm(current_strain)});
    else if(P == OutputType::PEEQ) data.emplace_back(vec{sqrt(2. / 3.) * tensor::strain::norm(current_strain - solve(initial_stiffness, current_stress))});
    else if(P == OutputType::SP) { if(vec principal_stress; eig_sym(principal_stress, tensor::stress::to_tensor(current_stress))) data.emplace_back(principal_stress); }
    else if(P == OutputType::EP) { if(vec principal_strain; eig_sym(principal_strain, tensor::strain::to_tensor(current_strain))) data.emplace_back(principal_strain); }
    else if(P == OutputType::EEP) { if(vec principal_strain; eig_sym(principal_strain, tensor::strain::to_tensor(solve(initial_stiffness, current_stress)))) data.emplace_back(principal_strain); }
    else if(P == OutputType::PEP) { if(vec principal_strain; eig_sym(principal_strain, tensor::strain::to_tensor(current_strain - solve(initial_stiffness, current_stress)))) data.emplace_back(principal_strain); }
    else if(P == OutputType::S11) data.emplace_back(vec{current_stress.n_elem > 0 ? current_stress(0) : 0.});
    else if(P == OutputType::S22) data.emplace_back(vec{current_stress.n_elem > 1 ? current_stress(1) : 0.});
    else if(P == OutputType::S33) data.emplace_back(vec{current_stress.n_elem > 2 ? current_stress(2) : 0.});
    else if(P == OutputType::S12) data.emplace_back(vec{current_stress.n_elem > 3 ? current_stress(3) : 0.});
    else if(P == OutputType::S23) data.emplace_back(vec{current_stress.n_elem > 4 ? current_stress(4) : 0.});
    else if(P == OutputType::S13) data.emplace_back(vec{current_stress.n_elem > 5 ? current_stress(5) : 0.});
    else if(P == OutputType::E11) data.emplace_back(vec{current_strain.n_elem > 0 ? current_strain(0) : 0.});
    else if(P == OutputType::E22) data.emplace_back(vec{current_strain.n_elem > 1 ? current_strain(1) : 0.});
    else if(P == OutputType::E33) data.emplace_back(vec{current_strain.n_elem > 2 ? current_strain(2) : 0.});
    else if(P == OutputType::E12) data.emplace_back(vec{current_strain.n_elem > 3 ? current_strain(3) : 0.});
    else if(P == OutputType::E23) data.emplace_back(vec{current_strain.n_elem > 4 ? current_strain(4) : 0.});
    else if(P == OutputType::E13) data.emplace_back(vec{current_strain.n_elem > 5 ? current_strain(5) : 0.});
    else if(P == OutputType::SP1) { if(vec principal_stress; eig_sym(principal_stress, tensor::stress::to_tensor(current_stress))) data.emplace_back(vec{principal_stress(0)}); }
    else if(P == OutputType::SP2) { if(vec principal_stress; eig_sym(principal_stress, tensor::stress::to_tensor(current_stress))) data.emplace_back(vec{principal_stress(1)}); }
    else if(P == OutputType::SP3) { if(vec principal_stress; eig_sym(principal_stress, tensor::stress::to_tensor(current_stress))) data.emplace_back(vec{principal_stress(2)}); }
    else if(P == OutputType::EP1) { if(vec principal_strain; eig_sym(principal_strain, tensor::strain::to_tensor(current_strain))) data.emplace_back(vec{principal_strain(0)}); }
    else if(P == OutputType::EP2) { if(vec principal_strain; eig_sym(principal_strain, tensor::strain::to_tensor(current_strain))) data.emplace_back(vec{principal_strain(1)}); }
    else if(P == OutputType::EP3) { if(vec principal_strain; eig_sym(principal_strain, tensor::strain::to_tensor(current_strain))) data.emplace_back(vec{principal_strain(2)}); }
    else if(P == OutputType::HIST) data.emplace_back(current_history);

    return data;
}
