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

#include "ElasticOS.h"

ElasticOS::ElasticOS(const unsigned T, const double E, const double P, const double R)
    : DataElasticOS{fabs(E), fabs(P)}
    , MaterialOS(T, R) {}

int ElasticOS::initialize(const shared_ptr<DomainBase>&) {
    const auto shear_modulus = elastic_modulus / (2. + 2. * poissons_ratio);
    initial_stiffness = diagmat(vec{elastic_modulus, shear_modulus, shear_modulus});

    ConstantStiffness(this);

    return SUANPAN_SUCCESS;
}

double ElasticOS::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return elastic_modulus;
    if(ParameterType::SHEARMODULUS == P || ParameterType::G == P) return elastic_modulus / (2. + 2. * poissons_ratio);
    if(ParameterType::BULKMODULUS == P) return elastic_modulus / (3. - 6. * poissons_ratio);
    if(ParameterType::POISSONSRATIO == P) return poissons_ratio;
    return 0.;
}

unique_ptr<Material> ElasticOS::get_copy() { return make_unique<ElasticOS>(*this); }

int ElasticOS::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;
    trial_stress = trial_stiffness * trial_strain;
    return SUANPAN_SUCCESS;
}

int ElasticOS::clear_status() {
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    return SUANPAN_SUCCESS;
}

int ElasticOS::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    return SUANPAN_SUCCESS;
}

int ElasticOS::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    return SUANPAN_SUCCESS;
}

void ElasticOS::print() {
    suanpan_info("An elastic material with E={:.4E} and nu={:.4E}.\n", elastic_modulus, poissons_ratio);
}
