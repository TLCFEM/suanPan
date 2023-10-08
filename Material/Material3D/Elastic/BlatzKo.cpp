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

#include "BlatzKo.h"
#include <Toolbox/tensor.h>

const vec BlatzKo::weight{2., 2., 2., 1., 1., 1.};

BlatzKo::BlatzKo(const unsigned T, const double E, const double V, const double R)
    : DataBlatzKo{fabs(E), fabs(V), fabs(E) / (2. + 2. * fabs(V)), (1. - fabs(V)) / (1. - 2. * fabs(V))}
    , Material3D(T, R) {}

int BlatzKo::initialize(const shared_ptr<DomainBase>&) {
    update_trial_status(zeros(6));
    current_stiffness = initial_stiffness = trial_stiffness;

    return SUANPAN_SUCCESS;
}

double BlatzKo::get_parameter(const ParameterType P) const { return material_property(elastic_modulus, poissons_ratio)(P); }

unique_ptr<Material> BlatzKo::get_copy() { return make_unique<BlatzKo>(*this); }

// takes green strain as input
int BlatzKo::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;

    vec G = weight % t_strain + tensor::unit_tensor2;

    vec H(6);
    H(0) = G(1) * G(2) - G(4) * G(4);
    H(1) = G(2) * G(0) - G(5) * G(5);
    H(2) = G(0) * G(1) - G(3) * G(3);
    H(3) = G(4) * G(5) - G(2) * G(3);
    H(4) = G(5) * G(3) - G(0) * G(4);
    H(5) = G(3) * G(4) - G(1) * G(5);

    const auto I3 = G(0) * H(0) + G(3) * H(3) + G(5) * H(5);

    auto factor_a = pow(std::max(datum::eps, I3), -half_beta_two);

    trial_stress = shear_modulus * (tensor::unit_tensor2 - factor_a * H);

    G *= -(factor_a *= 2. * shear_modulus);

    trial_stiffness.zeros(6, 6);

    trial_stiffness(4, 4) = -.5 * (trial_stiffness(1, 2) = G(0));
    trial_stiffness(5, 5) = -.5 * (trial_stiffness(0, 2) = G(1));
    trial_stiffness(3, 3) = -.5 * (trial_stiffness(0, 1) = G(2));
    trial_stiffness(4, 5) = -.5 * (trial_stiffness(2, 3) = -G(3));
    trial_stiffness(3, 5) = -.5 * (trial_stiffness(0, 4) = -G(4));
    trial_stiffness(3, 4) = -.5 * (trial_stiffness(1, 5) = -G(5));

    factor_a *= half_beta_two / I3;

    for(auto I = 0; I < 6; ++I) {
        const auto factor_b = factor_a * H(I);
        trial_stiffness(I, I) += factor_b * H(I);
        for(auto J = I + 1; J < 6; ++J) trial_stiffness(J, I) = trial_stiffness(I, J) += factor_b * H(J);
    }

    return SUANPAN_SUCCESS;
}

int BlatzKo::clear_status() {
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return SUANPAN_SUCCESS;
}

int BlatzKo::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int BlatzKo::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void BlatzKo::print() {
    suanpan_info("A Blatz-Ko material model.\n");
}
