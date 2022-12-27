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

#include "AxisymmetricElastic.h"
#include <Toolbox/tensorToolbox.h>

AxisymmetricElastic::AxisymmetricElastic(const unsigned T, const double E, const double P, const double R)
    : Material2D(T, PlaneType::A, R)
    , elastic_modulus(E)
    , poissons_ratio(P) {}

int AxisymmetricElastic::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio)(span(0, 3), span(0, 3));

    return SUANPAN_SUCCESS;
}

double AxisymmetricElastic::get_parameter(const ParameterType P) const {
    if(ParameterType::PLANETYPE == P) return static_cast<double>(plane_type);
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return elastic_modulus;
    if(ParameterType::SHEARMODULUS == P || ParameterType::G == P) return elastic_modulus / (2. + 2. * poissons_ratio);
    if(ParameterType::POISSONSRATIO == P) return poissons_ratio;
    return 0.;
}

unique_ptr<Material> AxisymmetricElastic::get_copy() { return make_unique<AxisymmetricElastic>(*this); }

int AxisymmetricElastic::update_trial_status(const vec& t_strain) {
    trial_stress = trial_stiffness * (trial_strain = t_strain);
    return SUANPAN_SUCCESS;
}

int AxisymmetricElastic::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    return reset_status();
}

int AxisymmetricElastic::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    return SUANPAN_SUCCESS;
}

int AxisymmetricElastic::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    return SUANPAN_SUCCESS;
}

void AxisymmetricElastic::print() {
    get_trial_strain().t().print("Strain:");
    get_trial_stress().t().print("Stress:");
}
