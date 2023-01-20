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

#include "Elastic1D.h"

Elastic1D::Elastic1D(const unsigned T, const double E, const double R)
    : DataElastic1D{E}
    , Material1D(T, R) {}

int Elastic1D::initialize(const shared_ptr<DomainBase>&) {
    initial_stiffness = elastic_modulus;
    ConstantStiffness(this);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Elastic1D::get_copy() { return make_unique<Elastic1D>(*this); }

int Elastic1D::update_trial_status(const vec& t_strain) {
    trial_stress = elastic_modulus * (trial_strain = t_strain);
    return SUANPAN_SUCCESS;
}

int Elastic1D::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    return 0;
}

int Elastic1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    return 0;
}

int Elastic1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    return 0;
}

void Elastic1D::print() {
    suanpan_info("A uniaxial elastic material with an elastic modulus of {:.4E}.\n", elastic_modulus);
    Material1D::print();
}
