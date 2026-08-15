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

#include "NonlocalIsotropicElastic3D.h"

#include <Toolbox/tensor.h>

const uvec NonlocalIsotropicElastic3D::u_dof{0, 1, 2, 3, 4, 5};
const uvec NonlocalIsotropicElastic3D::d_dof{6};

NonlocalIsotropicElastic3D::NonlocalIsotropicElastic3D(const unsigned T, const double E, const double P, const double ME, const double ER, const double R)
    : DataNonlocalIsotropicElastic3D{.elastic_modulus = std::fabs(E), .poissons_ratio = std::fabs(P), .maximum_energy = std::fabs(ME), .evolution_rate = std::fabs(ER / E)}
    , Material3D(T, R) {}

int NonlocalIsotropicElastic3D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio, nonlocal_size());

    initialize_history(1u);

    return SUANPAN_SUCCESS;
}

double NonlocalIsotropicElastic3D::get(const Parameter P) const { return MaterialProperty(elastic_modulus, poissons_ratio)(P); }

unique_ptr<Material> NonlocalIsotropicElastic3D::unique_copy() { return std::make_unique<NonlocalIsotropicElastic3D>(*this); }

unsigned NonlocalIsotropicElastic3D::nonlocal_size() const { return 1u; }

int NonlocalIsotropicElastic3D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;
    trial_stress = (trial_stiffness = initial_stiffness) * trial_strain;

    trial_history = current_history;

    if(const auto excessive_energy = .5 * dot(trial_stress, trial_strain) - maximum_energy; excessive_energy > trial_history(0)) {
        trial_history(0) = excessive_energy;
        trial_stress(6) = 1. - std::exp(-excessive_energy * evolution_rate);
        trial_stiffness(d_dof, u_dof) = (trial_stress(6) - 1.) * evolution_rate * trial_stress(u_dof).t();
    }
    else trial_stress(6) = current_stress(6);

    trial_stiffness(u_dof, d_dof) = trial_stress(u_dof);

    trial_stress(u_dof) *= 1. - trial_strain(6);
    trial_stiffness(u_dof, u_dof) *= 1. - trial_strain(6);

    return SUANPAN_SUCCESS;
}

int NonlocalIsotropicElastic3D::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    current_history = trial_history = initial_history;
    return SUANPAN_SUCCESS;
}

int NonlocalIsotropicElastic3D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    current_history = trial_history;
    return SUANPAN_SUCCESS;
}

int NonlocalIsotropicElastic3D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    trial_history = current_history;
    return SUANPAN_SUCCESS;
}

void NonlocalIsotropicElastic3D::print() {
    suanpan_info("A 3D isotropic elastic material with E={:.4E} and nu={:.4E}.\n", elastic_modulus, poissons_ratio);
}
