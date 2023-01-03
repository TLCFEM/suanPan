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

#include "IsotropicNonlinearElastic3D.h"
#include <Toolbox/tensorToolbox.h>

const mat IsotropicNonlinearElastic3D::unit_dev_tensor = two_third * tensor::unit_deviatoric_tensor4();
const mat IsotropicNonlinearElastic3D::unit_unit = unit_dev_tensor * tensor::unit_deviatoric_tensor4v2();

IsotropicNonlinearElastic3D::IsotropicNonlinearElastic3D(const unsigned T, const double R)
    : Material3D(T, R) {}

int IsotropicNonlinearElastic3D::initialize(const shared_ptr<DomainBase>&) {
    const auto derivative = compute_derivative(0., 0.);

    const auto& pwpd = derivative(1);
    const auto& ppwppm = derivative(2);

    const auto& pmpe = tensor::unit_tensor2;

    trial_stiffness = current_stiffness = initial_stiffness = ppwppm * pmpe * pmpe.t() + 2. * pwpd * unit_unit;

    return SUANPAN_SUCCESS;
}

int IsotropicNonlinearElastic3D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) < 1E-12) return SUANPAN_SUCCESS;

    const auto d_strain = tensor::dev(trial_strain);                                      // deviatoric strain
    const auto e_strain = two_third * dot(tensor::strain::norm_weight, square(d_strain)); // equivalent strain squared

    const auto derivative = compute_derivative(tensor::trace3(trial_strain), e_strain);

    const auto& pwpm = derivative(0);
    const auto& pwpd = derivative(1);
    const auto& ppwppm = derivative(2);
    const auto& ppwppd = derivative(3);
    const auto& ppwpmpd = derivative(4);
    const auto& ppwpdpm = derivative(5);

    const auto& pmpe = tensor::unit_tensor2;
    const vec pdpe = 2. * unit_dev_tensor * d_strain;

    trial_stress = pwpm * pmpe + pwpd * pdpe;

    trial_stiffness = pmpe * (ppwppm * pmpe.t() + ppwpmpd * pdpe.t()) + pdpe * (ppwppd * pdpe.t() + ppwpdpm * pmpe.t()) + 2. * pwpd * unit_unit;

    return SUANPAN_SUCCESS;
}

int IsotropicNonlinearElastic3D::clear_status() {
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return SUANPAN_SUCCESS;
}

int IsotropicNonlinearElastic3D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int IsotropicNonlinearElastic3D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void IsotropicNonlinearElastic3D::print() {
    suanpan_info("A 3D isotropic nonlinear elastic material model.\n");
}
