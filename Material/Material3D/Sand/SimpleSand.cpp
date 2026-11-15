/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
    : DataSimpleSand{.elastic_modulus = E, .poissons_ratio = V, .m = std::fabs(M), .a = A, .h = H, .ac = AC, .nb = std::fabs(NB), .nd = std::fabs(ND), .vc = std::fabs(VC), .pc = -std::fabs(PC), .lc = std::fabs(LC), .v0 = std::fabs(V0)}
    , Material3D(T, R) {}

int SimpleSand::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(6);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> SimpleSand::unique_copy() { return std::make_unique<SimpleSand>(*this); }

double SimpleSand::get(const Parameter P) const { return MaterialProperty(elastic_modulus, poissons_ratio)(P); }

int SimpleSand::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    const vec current_alpha(&current_history(0), 6);
    vec alpha(&trial_history(0), 6, false, true);

    const auto state_const = v0 - vc + v0 * tensor::trace<3>(trial_strain);

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    const auto trial_s = tensor::dev(trial_stress);
    const auto trial_p = tensor::mean<3>(trial_stress);

    // the yield surface and the critical state are defined for p<0 only, the sand has no
    // tensile strength; cap a non-negative trial pressure at zero and keep the purely
    // deviatoric elastic response
    if(trial_p > -std::fabs(pc) * tolerance) {
        trial_stress = trial_s;
        trial_stiffness = double_shear * unit_dev_tensor;

        return SUANPAN_SUCCESS;
    }

    auto s = trial_s;
    auto p = trial_p;

    // the critical state is undefined at p=0; bound the pressure away from zero so that the
    // state parameter and its derivatives remain finite if an iteration overshoots into tension,
    // the hard bound is replaced by a smooth soft minimum so that the local iteration does not
    // cycle around the kink when the pressure repeatedly crosses the bound
    const auto p_reg = std::min(-datum::eps, -1E-2 * std::fabs(pc));
    const auto p_width = .1 * std::fabs(p_reg);

    // if the local iterations stall because dilatancy drives the pressure towards tension, no
    // admissible plastic solution exists; the capped response is then used, the cap is only
    // allowed if the trial state or an iterate has reached the tension boundary
    auto at_bound = false;

    const auto tension_cap = [&] {
        if(!at_bound && trial_p <= -5E-2 * std::fabs(pc)) return false;

        trial_history = current_history;
        trial_stress = trial_s;
        trial_stiffness = double_shear * unit_dev_tensor;

        return true;
    };

    mat::fixed<14, 14> jacobian(fill::none);
    vec::fixed<14> residual(fill::none), incre(fill::none);

    jacobian(sa, sa) = 0.;

    vec n;
    double alpha_d, alpha_b;
    auto gamma = 0.;

    auto counter = 0u;

    // backtracking line search on the local residual: when a step increases the residual norm,
    // the previous state is restored and a halved step is taken instead, this damps the
    // oscillation that appears when the local system approaches the narrow tension wedge
    double relax{1.}, pre_norm{0.};
    auto pre_gamma = 0.;
    auto pre_p = 0.;
    vec pre_s, pre_alpha;

    const auto progress = [&](const double scale) {
        gamma -= scale * incre(sa);
        p -= scale * incre(sb);
        s -= scale * incre(sc);
        alpha -= scale * incre(sd);
    };

    while(true) {
        if(++counter > max_iteration) {
            if(tension_cap()) return SUANPAN_SUCCESS;

            return SUANPAN_FAIL;
        }

        const vec eta = s + p * alpha;
        const auto norm_eta = tensor::stress::norm(eta);
        n = eta / norm_eta;

        residual(sa) = norm_eta + m * p;

        if(1u == counter && residual(sa) < 0.) return SUANPAN_SUCCESS;

        // smooth soft minimum of p and p_reg, the pressure is bounded for the state parameter only
        const auto shifted = (p_reg - p) / p_width;
        const auto soft = shifted > 0. ? 1. / (1. + std::exp(-shifted)) : std::exp(shifted) / (1. + std::exp(shifted));
        const auto p_bounded = p_reg - p_width * (std::max(shifted, 0.) + std::log1p(std::exp(-std::fabs(shifted))));
        const auto dlp = soft * lc / p_bounded;
        const auto state = state_const + lc * std::log(p_bounded / pc);
        alpha_d = ac * std::exp(nd * state);
        alpha_b = ac * std::exp(-nb * state);
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
        jacobian(sb, sb) = 1. + bulk * a * gamma * (alpha_d * nd * dlp - (dot(alpha, unit_alpha) - alpha_n * alpha_n) / norm_eta);
        jacobian(sb, sc) = bulk * a * gamma / norm_eta * (alpha_n * unit_n.t() - unit_alpha.t());
        jacobian(sb, sd) = p * jacobian(sb, sc) - bulk * a * gamma * unit_n.t();

        jacobian(sc, sa) = double_shear * n;
        jacobian(sc, sb) = double_shear * gamma / norm_eta * (alpha - alpha_n * n);
        jacobian(sc, sc) = double_shear * gamma / norm_eta * eye(6, 6) - double_shear * gamma / norm_eta * n * unit_n.t();
        jacobian(sc, sd) = p * jacobian(sc, sc);
        jacobian(sc, sc) += eye(6, 6);

        jacobian(sd, sa) = h * alpha_b_m * n - h * alpha;
        jacobian(sd, sb) = gamma * h * alpha_b_m / norm_eta * alpha - gamma * h * (alpha_b_m * alpha_n / norm_eta + alpha_b * nb * dlp) * n;
        jacobian(sd, sc) = gamma * h * alpha_b_m / norm_eta * eye(6, 6) - gamma * h * alpha_b_m / norm_eta * n * unit_n.t();
        jacobian(sd, sd) = p * jacobian(sd, sc) - (gamma * h + 1.) * eye(6, 6);

        // if the previous step did not reduce the residual, restore the state and halve the step
        if(1u < counter && suanpan::inf_norm(residual) > pre_norm) {
            if(relax < 1E-4) {
                if(tension_cap()) return SUANPAN_SUCCESS;

                return SUANPAN_FAIL;
            }

            gamma = pre_gamma;
            p = pre_p;
            s = pre_s;
            alpha = pre_alpha;

            progress(relax *= .5);

            if(p > p_reg) at_bound = true;

            continue;
        }

        if(!solve(incre, jacobian, residual)) {
            if(tension_cap()) return SUANPAN_SUCCESS;

            return SUANPAN_FAIL;
        }

        const auto error = suanpan::inf_norm(incre);
        // both the increment and the residual are checked against the current stress scale so that
        // the tolerance remains meaningful when the stress itself is large
        const auto error_scale = 1. + std::max(std::fabs(p), suanpan::inf_norm(s));
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * error_scale || (counter > 5u && suanpan::inf_norm(residual) < tolerance * error_scale)) break;

        pre_norm = suanpan::inf_norm(residual);
        pre_gamma = gamma;
        pre_p = p;
        pre_s = s;
        pre_alpha = alpha;

        progress(relax = 1.);

        if(p > p_reg) at_bound = true;
    }

    trial_stress = s + p * tensor::unit_tensor2;

    mat::fixed<14, 6> left(fill::none), right(fill::none);

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

void SimpleSand::print() {
    suanpan_info("A simple sand model.\n");
}
