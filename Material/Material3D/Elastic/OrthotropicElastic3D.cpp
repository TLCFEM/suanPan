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

#include "OrthotropicElastic3D.h"
#include <Toolbox/tensorToolbox.h>

OrthotropicElastic3D::OrthotropicElastic3D(const unsigned T, vec&& E, vec&& P, const double R)
    : Material3D(T, R)
    , modulus(std::forward<vec>(E))
    , poissons_ratio(std::forward<vec>(P)) {}

int OrthotropicElastic3D::initialize(const shared_ptr<DomainBase>&) {
    initial_stiffness = tensor::orthotropic_stiffness(modulus, poissons_ratio);

    ConstantStiffness(this);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> OrthotropicElastic3D::get_copy() { return make_unique<OrthotropicElastic3D>(*this); }

int OrthotropicElastic3D::update_trial_status(const vec& t_strain) {
    trial_stress = trial_stiffness * (trial_strain = t_strain);
    return SUANPAN_SUCCESS;
}

int OrthotropicElastic3D::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    return SUANPAN_SUCCESS;
}

int OrthotropicElastic3D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    return SUANPAN_SUCCESS;
}

int OrthotropicElastic3D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    return SUANPAN_SUCCESS;
}

void OrthotropicElastic3D::print() { sp_info("A 3D orthotropic elastic material model.\n"); }
