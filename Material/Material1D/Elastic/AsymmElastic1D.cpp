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

#include "AsymmElastic1D.h"

AsymmElastic1D::AsymmElastic1D(const unsigned T, const double TE, const double CE, const double R)
    : DataAsymmElastic1D{TE, CE}
    , Material1D(T, R) {}

int AsymmElastic1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = t_elastic_modulus;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> AsymmElastic1D::get_copy() { return make_unique<AsymmElastic1D>(*this); }

int AsymmElastic1D::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;

    trial_stiffness = trial_strain(0) >= 0. ? t_elastic_modulus : c_elastic_modulus;

    trial_stress = trial_stiffness * trial_strain;

    return SUANPAN_SUCCESS;
}

int AsymmElastic1D::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return 0;
}

int AsymmElastic1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return 0;
}

int AsymmElastic1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return 0;
}

void AsymmElastic1D::print() {
    suanpan_info("A uniaxial elastic material with an elastic modulus of {:.4E} for tension and {:.4E} for compression.\n", t_elastic_modulus, c_elastic_modulus);
    Material1D::print();
}
