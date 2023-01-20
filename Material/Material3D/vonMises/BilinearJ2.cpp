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

#include "BilinearJ2.h"
#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

constexpr double BilinearJ2::two_third = 2. / 3.;
const double BilinearJ2::root_two_third = sqrt(two_third);
const mat BilinearJ2::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

BilinearJ2::BilinearJ2(const unsigned T, const double E, const double V, const double Y, const double H, const double B, const double R)
    : DataBilinearJ2{fabs(E), V, fabs(Y), H, B}
    , Material3D(T, R) {}

int BilinearJ2::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(7);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> BilinearJ2::get_copy() { return make_unique<BilinearJ2>(*this); }

double BilinearJ2::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return elastic_modulus;
    if(ParameterType::SHEARMODULUS == P || ParameterType::G == P) return shear_modulus;
    if(ParameterType::BULKMODULUS == P) return elastic_modulus / (3. - 6. * poissons_ratio);
    if(ParameterType::POISSONSRATIO == P) return poissons_ratio;
    return 0.;
}

int BilinearJ2::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= tolerance) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    auto& plastic_strain = trial_history(0);
    vec back_stress(&trial_history(1), 6, false, true);

    const vec shifted_stress = tensor::dev(trial_stress) - back_stress;

    const auto norm_shifted_stress = tensor::stress::norm(shifted_stress);

    const auto yield_surf = yield_stress + isotropic_modulus * plastic_strain;

    if(const auto yield_func = norm_shifted_stress - (yield_surf > 0. ? root_two_third * yield_surf : 0.); yield_func > 0.) {
        const auto tmp_a = double_shear + two_third * (kinematic_modulus + (yield_surf > 0. ? isotropic_modulus : 0.));
        const auto gamma = yield_func / tmp_a;

        auto tmp_b = double_shear * gamma / norm_shifted_stress;
        trial_stress -= tmp_b * shifted_stress;

        back_stress += two_third * kinematic_modulus * gamma / norm_shifted_stress * shifted_stress;
        plastic_strain += root_two_third * gamma;

        tmp_b *= double_shear;
        trial_stiffness += (tmp_b - square_double_shear / tmp_a) / norm_shifted_stress / norm_shifted_stress * shifted_stress * shifted_stress.t() - tmp_b * unit_dev_tensor;
    }

    return SUANPAN_SUCCESS;
}

int BilinearJ2::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int BilinearJ2::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int BilinearJ2::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

vector<vec> BilinearJ2::record(const OutputType P) {
    if(P == OutputType::MISES) return {vec{tensor::stress::norm(tensor::dev(current_stress)) / root_two_third}};
    if(P == OutputType::EEQ) return {vec{root_two_third * tensor::strain::norm(tensor::dev(current_strain))}};
    if(P == OutputType::PEEQ) return {vec{current_history(0)}};

    return Material3D::record(P);
}

void BilinearJ2::print() {
    suanpan_info("A 3D bilinear hardening model.\n");
}
