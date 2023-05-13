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

#include "VAFCRP.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

const double VAFCRP::root_three_two = sqrt(1.5);
const mat VAFCRP::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

VAFCRP::VAFCRP(const unsigned T, const double E, const double V, const double Y, const double S, const double H, const double M, const double MU, const double EP, vec&& A, vec&& B, const double R)
    : DataVAFCRP{E, V, Y, S, H, M, MU, EP, std::forward<vec>(A), std::forward<vec>(B)}
    , Material3D(T, R) { access::rw(tolerance) = 1E-15; }

int VAFCRP::initialize(const shared_ptr<DomainBase>& D) {
    incre_time = &D->get_factory()->modify_incre_time();

    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(1 + 6 * size);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> VAFCRP::get_copy() { return make_unique<VAFCRP>(*this); }

double VAFCRP::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return elastic_modulus;
    if(ParameterType::SHEARMODULUS == P || ParameterType::G == P) return shear;
    if(ParameterType::BULKMODULUS == P) return elastic_modulus / (3. - 6. * poissons_ratio);
    if(ParameterType::POISSONSRATIO == P) return poissons_ratio;
    return 0.;
}

int VAFCRP::update_trial_status(const vec& t_strain) {
    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * (incre_strain = (trial_strain = t_strain) - current_strain);

    trial_history = current_history;
    auto& p = trial_history(0);

    const auto trial_s = tensor::dev(trial_stress);

    auto eta = trial_s;
    for(unsigned I = 0; I < size; ++I) eta -= vec{&trial_history(1 + 6llu * I), 6, false, true};

    // const auto residual = root_three_two * tensor::stress::norm(eta) - std::max(0., yield + hardening * p + saturated * (1. - exp(-m * p)));

    if(root_three_two * tensor::stress::norm(eta) < std::max(0., yield + hardening * p + saturated * (1. - exp(-m * p)))) return SUANPAN_SUCCESS;

    vec xi;
    auto gamma = 0., exp_gamma = 1.;
    double norm_xi, jacobian;

    unsigned counter = 0;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto exp_term = saturated * exp(-m * p);

        auto k = yield + saturated + hardening * p - exp_term;
        auto dk = 0.;
        k < 0. ? k = 0. : dk = hardening + m * exp_term;

        vec sum_a(6, fill::zeros);
        auto sum_b = 0.;
        for(unsigned I = 0; I < size; ++I) {
            const auto denom = 1. + b(I) * gamma;
            sum_a += vec{&trial_history(1 + 6llu * I), 6, false, true} / denom;
            sum_b += a(I) * gamma / denom;
        }

        norm_xi = tensor::stress::norm(xi = trial_s - sum_a);

        const auto q = root_three_two * (norm_xi - root_six_shear * gamma - sum_b);
        exp_gamma = pow(*incre_time / (*incre_time + mu * gamma), epsilon);

        sum_b = 0.;
        for(unsigned I = 0; I < size; ++I) sum_b += (b(I) / norm_xi * tensor::stress::double_contraction(xi, vec{&trial_history(1 + 6llu * I), 6, false, true}) - a(I)) * pow(1. + b(I) * gamma, -2.);

        jacobian = exp_gamma * (root_three_two * sum_b - three_shear - q * epsilon * mu / (*incre_time + mu * gamma)) - dk;

        const auto incre = (q * exp_gamma - k) / jacobian;
        suanpan_debug("Local iteration error: {:.5E}.\n", fabs(incre));
        if(fabs(incre) <= tolerance) break;

        gamma -= incre;
        p -= incre;
    }

    const vec u = xi / norm_xi;

    vec sum_c(6, fill::zeros);
    for(unsigned I = 0; I < size; ++I) {
        vec beta(&trial_history(1 + 6llu * I), 6, false, true);
        sum_c += b(I) * pow(1. + b(I) * gamma, -2.) * (beta - tensor::stress::double_contraction(u, beta) * u);
        beta = (beta + a(I) * gamma * u) / (1. + b(I) * gamma);
    }

    trial_stress -= root_six_shear * gamma * u;

    trial_stiffness += (root_six_shear * (double_shear * gamma / norm_xi + root_six_shear * exp_gamma / jacobian) * u + root_six_shear * root_six_shear * exp_gamma * gamma / jacobian / norm_xi * sum_c) * u.t() - double_shear * root_six_shear * gamma / norm_xi * unit_dev_tensor;

    return SUANPAN_SUCCESS;
}

int VAFCRP::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int VAFCRP::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int VAFCRP::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

vector<vec> VAFCRP::record(const OutputType P) {
    if(P == OutputType::PEEQ) return {vec{current_history(0)}};

    return Material3D::record(P);
}

void VAFCRP::print() {
    suanpan_info("A VADCRP material model.\n");
}
