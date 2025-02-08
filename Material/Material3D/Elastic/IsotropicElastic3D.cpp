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

#include "IsotropicElastic3D.h"
#include <Toolbox/tensor.h>

IsotropicElastic3D::IsotropicElastic3D(const unsigned T, const double E, const double P, const double R)
    : DataIsotropicElastic3D{fabs(E), fabs(P)}
    , Material3D(T, R) {}

int IsotropicElastic3D::initialize(const shared_ptr<DomainBase>&) {
    initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    ConstantStiffness(this);

    return SUANPAN_SUCCESS;
}

double IsotropicElastic3D::get_parameter(const ParameterType P) const { return material_property(elastic_modulus, poissons_ratio)(P); }

unique_ptr<Material> IsotropicElastic3D::get_copy() { return make_unique<IsotropicElastic3D>(*this); }

int IsotropicElastic3D::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;
    trial_stress = trial_stiffness * trial_strain;
    // incre_strain = trial_strain - current_strain;
    // incre_stress = trial_stress - current_stress;
    return SUANPAN_SUCCESS;
}

int IsotropicElastic3D::clear_status() {
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    return SUANPAN_SUCCESS;
}

int IsotropicElastic3D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    return SUANPAN_SUCCESS;
}

int IsotropicElastic3D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    return SUANPAN_SUCCESS;
}

void IsotropicElastic3D::print() {
    suanpan_info("A 3D isotropic elastic material with E={:.4E} and nu={:.4E}.\n", elastic_modulus, poissons_ratio);
}
