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

#include "PhaseField.h"
#include <Toolbox/tensor.h>

void PhaseField::commit_status(const unique_ptr<Material>& mat_obj) {
    const auto i_strain = mat_obj->get_trial_strain() - mat_obj->get_current_strain();

    if(norm(i_strain) <= 1E-14) return;

    vec D;
    mat P;

    if(const vec a_stress = .5 * (mat_obj->get_trial_stress() + mat_obj->get_current_stress()); eig_sym(D, P, tensor::stress::to_tensor(a_stress))) {
        vec PD = zeros(D.n_elem), ND = zeros(D.n_elem);
        for(auto I = 0llu; I < D.n_elem; ++I) (D(I) >= 0. ? PD(I) : ND(I)) = D(I);

        const auto t_stress = tensor::stress::to_voigt(P * diagmat(PD) * P.t());
        const auto c_stress = tensor::stress::to_voigt(P * diagmat(ND) * P.t());

        const auto t_dev_stress = tensor::dev(t_stress);
        const auto t_vol_stress = t_stress - t_dev_stress;
        const auto c_dev_stress = tensor::dev(c_stress);
        const auto c_vol_stress = c_stress - c_dev_stress;

        strain_energy += dot(i_strain, a_stress);
        t_strain_energy += dot(i_strain, t_stress);
        c_strain_energy += dot(i_strain, c_stress);
        dev_strain_energy += dot(i_strain, t_dev_stress + c_dev_stress);
        t_dev_strain_energy += dot(i_strain, t_dev_stress);
        c_dev_strain_energy += dot(i_strain, c_dev_stress);
        vol_strain_energy += dot(i_strain, t_vol_stress + c_vol_stress);
        t_vol_strain_energy += dot(i_strain, t_vol_stress);
        c_vol_strain_energy += dot(i_strain, c_vol_stress);
    }
}

void PhaseField::clear_status() {
    strain_energy = 0.;
    t_strain_energy = 0.;
    c_strain_energy = 0.;
    dev_strain_energy = 0.;
    t_dev_strain_energy = 0.;
    c_dev_strain_energy = 0.;
    vol_strain_energy = 0.;
    t_vol_strain_energy = 0.;
    c_vol_strain_energy = 0.;
}
