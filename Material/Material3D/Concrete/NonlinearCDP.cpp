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

#include "NonlinearCDP.h"

#include <Recorder/OutputType.h>
#include <Toolbox/ridders.hpp>
#include <Toolbox/tensor.h>

const double NonlinearCDP::root_three_two = std::sqrt(1.5);
const mat NonlinearCDP::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

double NonlinearCDP::compute_r(const vec3& in) {
    const auto r = .5 + .5 * accu(in) / accu(abs(in));
    return !std::isfinite(r) || r < 0. ? .0 : r > 1. ? 1. : r;
}

vec3 NonlinearCDP::compute_dr(const vec3& in) {
    const auto g = accu(abs(in));

    vec3 out = .5 / g * (ones(3) - accu(in) * sign(in) / g);

    if(!out.is_finite()) out.zeros();

    return out;
}

double NonlinearCDP::compute_s(const double r) const { return s0 + r - s0 * r; }

NonlinearCDP::NonlinearCDP(const unsigned T, const double E, const double V, const double GT, const double GC, const double AP, const double BC, const double S, const double R)
    : DataNonlinearCDP{std::fabs(E), V < .5 ? V : .2, std::fabs(GT), std::fabs(GC), (std::fabs(BC) - 1.) / (2. * std::fabs(BC) - 1.), std::fabs(AP), std::fabs(S)}
    , Material3D(T, R) { access::rw(tolerance) = 1E-13; }

int NonlinearCDP::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(10);

    return SUANPAN_SUCCESS;
}

double NonlinearCDP::get_parameter(const ParameterType P) const { return material_property(elastic_modulus, poissons_ratio)(P); }

int NonlinearCDP::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& d_t = trial_history(0);
    auto& d_c = trial_history(1);
    auto& kappa_t = trial_history(2);
    auto& kappa_c = trial_history(3);
    vec plastic_strain(&trial_history(4), 6, false, true);

    const auto& current_kappa_t = current_history(2);
    const auto& current_kappa_c = current_history(3);

    const auto bound_kappa_t = [&] {
        if(kappa_t > 1.) kappa_t = 1. - std::numeric_limits<float>::epsilon(); // avoid overshoot
        else if(kappa_t < current_kappa_t) kappa_t = current_kappa_t;
        return kappa_t;
    };
    const auto bound_kappa_c = [&] {
        if(kappa_c > 1.) kappa_c = 1. - std::numeric_limits<float>::epsilon(); // avoid overshoot
        else if(kappa_c < current_kappa_c) kappa_c = current_kappa_c;
        return kappa_c;
    };

    trial_stress = (trial_stiffness = initial_stiffness) * (trial_strain - plastic_strain); // 6

    vec3 principal_stress;     // 3
    mat33 principal_direction; // 3x3
    if(!eig_sym(principal_stress, principal_direction, tensor::stress::to_tensor(trial_stress), "std")) return SUANPAN_FAIL;

    const auto trans = transform::compute_jacobian_nominal_to_principal(principal_direction);

    const auto s = tensor::dev(trial_stress);    // 6
    const auto norm_s = tensor::stress::norm(s); // 1
    vec6 n = s / norm_s;                         // 6
    if(!n.is_finite()) n.zeros();

    const auto ps = tensor::dev(principal_stress); // 3
    const vec3 pn = normalise(ps);                 // 3

    const vec3 dsigmadlambda = -double_shear * pn - three_alpha_p_bulk; // 3

    const auto dgdsigma_t = (pn(2) + alpha_p) / g_t;
    const auto dgdsigma_c = (pn(0) + alpha_p) / g_c;

    auto new_stress = principal_stress; // converged principal stress

    const auto const_yield = alpha * accu(principal_stress) + root_three_two * norm_s;

    vec3 residual, incre;
    mat33 jacobian(fill::zeros);

    pod6 t_para, c_para;

    auto lambda = 0., ref_error = 0.;
    double r, beta;
    vec3 dr;

    auto counter = 0u;
    auto try_bisection = false;
    while(true) {
        if(max_iteration == ++counter) {
            if(try_bisection) {
                suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
                return SUANPAN_FAIL;
            }

            try_bisection = true;
            counter = 2u; // skip elastic check

            auto approx_kappa_t = [&](const double in_kappa) {
                t_para = compute_tension_backbone(kappa_t = in_kappa);
                return r * t_para[1] * dgdsigma_t * lambda + current_kappa_t - kappa_t;
            };
            auto approx_kappa_c = [&](const double in_kappa) {
                c_para = compute_compression_backbone(kappa_c = in_kappa);
                return (1. - r) * c_para[1] * dgdsigma_c * lambda + current_kappa_c - kappa_c;
            };

            auto approx_update = [&](const double in_lambda) {
                r = compute_r(new_stress = principal_stress + (lambda = in_lambda) * dsigmadlambda);

                ridders(approx_kappa_t, current_kappa_t, 1. - datum::eps, tolerance);
                ridders(approx_kappa_c, current_kappa_c, 1. - datum::eps, tolerance);

                auto f = const_yield + pfplambda * lambda + one_minus_alpha * c_para[2];

                if(new_stress(2) > 0.) f -= (one_minus_alpha * c_para[2] / t_para[2] + alpha + 1.) * new_stress(2);

                return f;
            };

            auto x1{0.};
            lambda = .25 * tensor::strain::norm(incre_strain) / std::sqrt(1. + 3. * alpha_p * alpha_p);
            while(approx_update(lambda) > 0.) {
                x1 = lambda;
                lambda *= 2.;
            }

            ridders(approx_update, x1, lambda, tolerance);
        }

        t_para = compute_tension_backbone(kappa_t);
        c_para = compute_compression_backbone(kappa_c);

        beta = -one_minus_alpha * c_para[2] / t_para[2] - alpha - 1.;

        residual(0) = const_yield + pfplambda * lambda + one_minus_alpha * c_para[2];

        if(new_stress(2) > 0.) residual(0) += beta * new_stress(2);

        r = compute_r(new_stress);

        if(1u == counter && residual(0) < 0.) {
            const auto damage_c = scale * d_c - 1.;
            const auto damage_t = compute_s(r) * scale * d_t - 1.;
            const auto damage = damage_c * damage_t;
            trial_stiffness = (damage * eye(6, 6) + damage_c * scale * d_t * (1. - s0) * trial_stress * compute_dr(new_stress).t() * trans) * initial_stiffness;
            trial_stress *= damage;
            return SUANPAN_SUCCESS;
        }

        const auto t_term = t_para[1] * dgdsigma_t;
        const auto c_term = c_para[1] * dgdsigma_c;

        residual(1) = r * t_term * lambda + current_kappa_t - kappa_t;
        residual(2) = (c_term - r * c_term) * lambda + current_kappa_c - kappa_c;

        if(new_stress(2) > 0.) {
            jacobian(0, 0) = pfplambda + beta * dsigmadlambda(2);
            const auto tmp_term = one_minus_alpha * new_stress(2) / t_para[2];
            jacobian(0, 1) = tmp_term * c_para[2] / t_para[2] * t_para[5];
            jacobian(0, 2) = (one_minus_alpha - tmp_term) * c_para[5];
        }
        else {
            jacobian(0, 0) = pfplambda;
            jacobian(0, 1) = 0.;
            jacobian(0, 2) = one_minus_alpha * c_para[5];
        }

        const auto dlambda = r + lambda * dot(dr = compute_dr(new_stress), dsigmadlambda);
        jacobian(1, 0) = t_term * dlambda;
        jacobian(2, 0) = c_term - c_term * dlambda;
        jacobian(1, 1) = r * lambda * dgdsigma_t * t_para[4] - 1.;
        jacobian(2, 2) = (lambda - r * lambda) * dgdsigma_c * c_para[4] - 1.;

        if(rcond(jacobian) < datum::eps) {
            // short-circuit and direct go to the bisection logic
            counter = max_iteration - 1u;
            continue;
        }

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) break;

        lambda -= incre(0);
        kappa_t -= incre(1);
        kappa_c -= incre(2);
        new_stress -= dsigmadlambda * incre(0);

        bound_kappa_t();
        bound_kappa_c();
    }

    // update damage indices
    d_t = t_para[0];
    d_c = c_para[0];
    // update plastic strain
    plastic_strain += lambda * (n % tensor::stress::norm_weight + unit_alpha_p);

    const auto recovery = compute_s(r);
    const auto damage_c = scale * d_c - 1.;
    const auto damage_t = recovery * scale * d_t - 1.;
    const auto damage = damage_c * damage_t;

    // update trial stress
    trial_stress = transform::compute_jacobian_principal_to_nominal(principal_direction) * new_stress;

    const mat66 dnde = double_shear / norm_s * (unit_dev_tensor - n * n.t());

    // \dfrac{\partial\bar{\sigma}}{\partial\varepsilon^{tr}}
    trial_stiffness -= double_shear * lambda * dnde;

    const rowvec6 drdsigma = dr.t() * trans;
    const rowvec6 prpe = drdsigma * trial_stiffness;

    // compute local derivatives
    mat left(3, 6);
    left.row(0) = 3. * alpha * bulk * tensor::unit_tensor2.t() + root_three_two * double_shear * n.t();
    left.row(1) = t_para[1] * lambda * (r / g_t * trans.row(2) * dnde + dgdsigma_t * prpe);
    left.row(2) = c_para[1] * lambda * ((1. - r) / g_c * trans.row(0) * dnde - dgdsigma_c * prpe);

    if(new_stress(2) > 0.) left.row(0) += beta * trans.row(2) * trial_stiffness;

    const mat right = -solve(jacobian, left);
    const auto& dlambdade = right.row(0);
    const auto& dkappade = right.rows(1, 2);

    // \dfrac{\mathrm{d}\bar{\sigma}}{\mathrm{d}\varepsilon^{tr}}
    trial_stiffness -= (double_shear * n + three_alpha_p_bulk * tensor::unit_tensor2) * dlambdade;

    trial_stiffness = (damage * eye(6, 6) + scale * d_t * damage_c * (1. - s0) * trial_stress * drdsigma) * trial_stiffness + trial_stress * scale * rowvec2{recovery * damage_c * t_para[3], damage_t * c_para[3]} * dkappade;

    trial_stress *= damage;

    return SUANPAN_SUCCESS;
}

int NonlinearCDP::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int NonlinearCDP::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int NonlinearCDP::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

std::vector<vec> NonlinearCDP::record(const OutputType P) {
    if(P == OutputType::DT) return {vec{current_history(0)}};
    if(P == OutputType::DC) return {vec{current_history(1)}};

    return Material3D::record(P);
}

void NonlinearCDP::print() {
    suanpan_info("A concrete damage plasticity model.\n");
}
