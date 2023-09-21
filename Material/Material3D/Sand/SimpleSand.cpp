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

#include "SimpleSand.h"
#include <Toolbox/tensor.h>

const span SimpleSand::sc(2, 7);
const span SimpleSand::sd(8, 13);
const mat SimpleSand::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

SimpleSand::SimpleSand(const unsigned T, const double E, const double V, const double M, const double A, const double H, const double AC, const double NB, const double ND, const double VC, const double PC, const double LC, const double V0, const double R)
    : DataSimpleSand{E, V, fabs(M), A, H, AC, fabs(NB), fabs(ND), fabs(VC), -fabs(PC), fabs(LC), fabs(V0)}
    , Material3D(T, R) { access::rw(tolerance) = 1E-12; }

int SimpleSand::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(6);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> SimpleSand::get_copy() { return make_unique<SimpleSand>(*this); }

double SimpleSand::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return elastic_modulus;
    if(ParameterType::SHEARMODULUS == P || ParameterType::G == P) return shear;
    if(ParameterType::BULKMODULUS == P) return elastic_modulus / (3. - 6. * poissons_ratio);
    if(ParameterType::POISSONSRATIO == P) return poissons_ratio;
    return 0.;
}

int SimpleSand::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    const vec current_alpha(&current_history(0), 6);
    vec alpha(&trial_history(0), 6, false, true);

    const auto state_const = v0 - vc + v0 * tensor::trace3(trial_strain);

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    const auto trial_s = tensor::dev(trial_stress);
    const auto trial_p = tensor::mean3(trial_stress);
    auto s = trial_s;
    auto p = trial_p;

    mat jacobian(14, 14, fill::none);
    vec residual(14, fill::none), incre;

    jacobian(sa, sa) = 0.;

    vec n;
    double alpha_d, alpha_b;
    auto gamma = 0., ref_error = 1.;

    auto counter = 0u;

    while(true) {
        if(max_iteration == ++counter) return SUANPAN_FAIL;

        const vec eta = s + p * alpha;
        const auto norm_eta = tensor::stress::norm(eta);
        n = eta / norm_eta;

        residual(sa) = norm_eta + m * p;

        if(1u == counter && residual(sa) < 0.) return SUANPAN_SUCCESS;

        const auto state = state_const + lc * log(p / pc);
        alpha_d = ac * exp(nd * state);
        alpha_b = ac * exp(-nb * state);
        const auto alpha_d_m = alpha_d - m;
        const auto alpha_b_m = alpha_b - m;
        const vec unit_n = n % tensor::stress::norm_weight;
        const vec unit_alpha = alpha % tensor::stress::norm_weight;
        const auto alpha_n = dot(unit_alpha, n);

        residual(sb) = p - trial_p + bulk * a * gamma * (alpha_d_m - alpha_n);
        residual(sc) = s - trial_s + double_shear * gamma * n;
        residual(sd) = current_alpha + gamma * h * alpha_b_m * n - (gamma * h + 1.) * alpha;

        jacobian(sa, sb) = m + alpha_n;
        jacobian(sa, sc) = unit_n.t();
        jacobian(sa, sd) = p * jacobian(sa, sc);

        jacobian(sb, sa) = bulk * a * (alpha_d_m - alpha_n);
        jacobian(sb, sb) = 1. + bulk * a * gamma * (alpha_d * nd * lc / p - (dot(alpha, unit_alpha) - alpha_n * alpha_n) / norm_eta);
        jacobian(sb, sc) = bulk * a * gamma / norm_eta * (alpha_n * unit_n.t() - unit_alpha.t());
        jacobian(sb, sd) = p * jacobian(sb, sc) - bulk * a * gamma * unit_n.t();

        jacobian(sc, sa) = double_shear * n;
        jacobian(sc, sb) = double_shear * gamma / norm_eta * (alpha - alpha_n * n);
        jacobian(sc, sc) = double_shear * gamma / norm_eta * eye(6, 6) - double_shear * gamma / norm_eta * n * unit_n.t();
        jacobian(sc, sd) = p * jacobian(sc, sc);
        jacobian(sc, sc) += eye(6, 6);

        jacobian(sd, sa) = h * alpha_b_m * n - h * alpha;
        jacobian(sd, sb) = gamma * h * alpha_b_m / norm_eta * alpha - gamma * h * (alpha_b_m * alpha_n / norm_eta + alpha_b * nb * lc / p) * n;
        jacobian(sd, sc) = gamma * h * alpha_b_m / norm_eta * eye(6, 6) - gamma * h * alpha_b_m / norm_eta * n * unit_n.t();
        jacobian(sd, sd) = p * jacobian(sd, sc) - (gamma * h + 1.) * eye(6, 6);

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        auto error = norm(residual);

        if(1u == counter) ref_error = std::max(1., error);
        suanpan_debug("Local iteration error: {:.5E}.\n", error /= ref_error);
        if(error <= tolerance || norm(incre) <= tolerance) break;

        gamma -= incre(sa);
        p -= incre(sb);
        s -= incre(sc);
        alpha -= incre(sd);
    }

    trial_stress = s + p * tensor::unit_tensor2;

    mat left(14, 6, fill::none), right;

    left.row(sa).zeros();
    left.row(sb) = (bulk - bulk * a * gamma * alpha_d * nd * v0) * tensor::unit_tensor2.t();
    left.rows(sc) = double_shear * unit_dev_tensor;
    left.rows(sd) = alpha_b * nb * v0 * gamma * h * n * tensor::unit_tensor2.t();

    if(!solve(right, jacobian, left)) return SUANPAN_FAIL;

    trial_stiffness = right.rows(sc);
    trial_stiffness.row(0) += right.row(sb);
    trial_stiffness.row(1) += right.row(sb);
    trial_stiffness.row(2) += right.row(sb);

    return SUANPAN_SUCCESS;
}

int SimpleSand::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int SimpleSand::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int SimpleSand::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void SimpleSand::print() {
    suanpan_info("A simple sand model.\n");
}
