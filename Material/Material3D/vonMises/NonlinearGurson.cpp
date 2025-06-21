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

#include "NonlinearGurson.h"

#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

const double NonlinearGurson::sqrt_three_two = std::sqrt(1.5);
const mat NonlinearGurson::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

NonlinearGurson::NonlinearGurson(const unsigned T, const double E, const double V, const double Q1, const double Q2, const double FN, const double SN, const double EN, const double R)
    : DataNonlinearGurson{E, V, Q1, Q2, FN, SN, EN}
    , Material3D(T, R) {}

int NonlinearGurson::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(2);

    return SUANPAN_SUCCESS;
}

double NonlinearGurson::get_parameter(const ParameterType P) const { return material_property(elastic_modulus, poissons_ratio)(P); }

int NonlinearGurson::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    auto& pe = trial_history(0); // equivalent plastic strain
    auto& f = trial_history(1);  // volume fraction
    const auto& current_pe = current_history(0);
    const auto& current_f = current_history(1);

    auto trial_s = tensor::dev(trial_stress);                            // trial deviatoric stress
    const auto trial_q = sqrt_three_two * tensor::stress::norm(trial_s); // trial von Mises stress
    const auto trial_p = tensor::mean3(trial_stress);                    // trial hydrostatic stress
    auto p = trial_p;                                                    // hydrostatic stress

    mat44 jacobian(fill::none);
    vec4 incre(fill::none), residual(fill::none);
    auto gamma = 0.;
    double denom;

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto hardening = compute_hardening(pe);
        const auto &k = hardening(0), &dk = hardening(1);
        const auto hyper_term = 1.5 * q2 * p / k;
        const auto cosh_term = std::cosh(hyper_term);
        const auto sinh_term = std::sinh(hyper_term);
        const auto q = trial_q / (denom = 1. + six_shear * gamma);
        const auto an = para_b * std::exp(-.5 * std::pow((pe - en) / sn, 2.));
        const auto para_d = para_a * sinh_term;

        const auto diff_pe = pe - current_pe, diff_p = p - trial_p;

        residual(0) = q * q + k * k * (f * q1 * (2. * cosh_term - q1 * f) - 1.);

        if(1 == counter && residual(0) < 0.) return SUANPAN_SUCCESS;

        residual(1) = (1. - f) * k * diff_pe - 2. * gamma * q * q + p * diff_p / bulk;
        residual(2) = f - current_f + (1. - f) * diff_p / bulk - an * diff_pe;
        residual(3) = diff_p + para_a * gamma * f * k * sinh_term;

        jacobian(0, 0) = -2. * six_shear / denom * q * q;
        jacobian(0, 1) = (f * (4. * q1 * k * cosh_term - para_d / bulk * p) - 2. * k * (q1 * q1 * f * f + 1.)) * dk;
        jacobian(0, 2) = 2. * k * k * q1 * (cosh_term - q1 * f);
        jacobian(0, 3) = para_d / bulk * f * k;
        jacobian(1, 0) = 2. * q * q * (six_shear * gamma - 1.) / denom;
        jacobian(1, 1) = (1. - f) * (dk * diff_pe + k);
        jacobian(1, 2) = -k * diff_pe;
        jacobian(1, 3) = (p + diff_p) / bulk;
        jacobian(2, 0) = 0.;
        jacobian(2, 1) = an / sn / sn * (pe - en) * diff_pe - an;
        jacobian(2, 2) = 1. - diff_p / bulk;
        jacobian(2, 3) = (1. - f) / bulk;
        jacobian(3, 0) = para_d * f * k;
        jacobian(3, 1) = para_a * gamma * f * (sinh_term - hyper_term * cosh_term) * dk;
        jacobian(3, 2) = para_d * gamma * k;
        jacobian(3, 3) = 1. + 1.5 * para_a * q2 * gamma * f * cosh_term;

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) break;

        gamma -= incre(0);
        pe -= incre(1);
        f -= incre(2);
        p -= incre(3);

        f = std::min(std::max(f, 0.), 1.); // avoid overshoot
        pe = std::max(pe, 0.);
    }

    trial_s /= denom;

    mat::fixed<4, 6> left(fill::none), right(fill::none);

    right.row(0) = -six_shear / denom * trial_s.t();
    right.row(1) = -2. * gamma * right.row(0) + p * tensor::unit_tensor2.t();
    right.row(2) = (1. - f) * tensor::unit_tensor2.t();
    right.row(3) = bulk * tensor::unit_tensor2.t();

    if(!solve(left, jacobian, right)) return SUANPAN_FAIL;

    trial_stress = trial_s + p * tensor::unit_tensor2;

    trial_stiffness = six_shear / denom / 3. * unit_dev_tensor - six_shear / denom * trial_s * left.row(0) + tensor::unit_tensor2 * left.row(3);

    return SUANPAN_SUCCESS;
}

int NonlinearGurson::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int NonlinearGurson::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int NonlinearGurson::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

std::vector<vec> NonlinearGurson::record(const OutputType P) {
    if(P == OutputType::VF) return {vec{current_history(1)}};

    return Material3D::record(P);
}

void NonlinearGurson::print() {
    suanpan_info("A 3D nonlinear hardening model using Gurson's model.\n");
}
