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

#include "Subloading.h"

#include <Toolbox/tensor.h>

const double Subloading::root_two_third = sqrt(two_third);
const double Subloading::rate_bound = -log(z_bound);
const mat Subloading::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

pod2 Subloading::yield_ratio(const double z) {
    if(z < z_bound) return {rate_bound, 0.};

    return {-log(z), -1. / z};
}

Subloading::Subloading(const unsigned T, DataSubloading&& D, const double R)
    : DataSubloading{std::move(D)}
    , Material3D(T, R) { access::rw(tolerance) = 1E-13; }

int Subloading::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic, poissons_ratio);

    initialize_history(15);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Subloading::get_copy() { return make_unique<Subloading>(*this); }

int Subloading::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    const auto& current_q = current_history(1);
    const auto& current_z = current_history(2);
    const vec current_alpha(&current_history(3), 6, false, true);
    const vec current_d(&current_history(9), 6, false, true);
    auto& iteration = trial_history(0);
    auto& q = trial_history(1);
    auto& z = trial_history(2);
    vec alpha(&trial_history(3), 6, false, true);
    vec d(&trial_history(9), 6, false, true);

    const vec trial_s = tensor::dev(trial_stress);

    iteration = 0.;
    auto gamma = 0., ref_error = 0.;
    auto start_z = current_z;

    vec2 residual, incre;
    mat22 jacobian;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        q = current_q + root_two_third * gamma;

        const auto exp_iso = saturation_iso * exp(-m_iso * q);
        auto y = initial_iso + saturation_iso + k_iso * q - exp_iso;
        auto dy = root_two_third * (k_iso + m_iso * exp_iso);
        if(y < 0.) y = dy = 0.;

        const auto exp_kin = saturation_kin * exp(-m_kin * q);
        auto a = initial_kin + saturation_kin + k_kin * q - exp_kin;
        auto da = root_two_third * (k_kin + m_kin * exp_kin);
        if(a < 0.) a = da = 0.;

        const auto bot_alpha = 1. + b.r() * gamma;
        const auto bot_d = 1. + c.r() * gamma;

        const vec pzetapz = y / bot_d * current_d;
        const vec pzetapgamma = (b.r() * a / bot_alpha - da) / bot_alpha * current_alpha + (z - 1.) * (dy - c.r() * y / bot_d) / bot_d * current_d;

        const vec zeta = trial_s - a / bot_alpha * current_alpha + (z - 1.) * pzetapz;
        const auto norm_zeta = tensor::stress::norm(zeta);
        const vec n = zeta / norm_zeta;

        alpha = (root_two_third * gamma * b.rb() * n + current_alpha) / bot_alpha;

        d = (root_two_third * gamma * c.rb() * n + current_d) / bot_d;

        if(1u == counter) {
            const vec ref = trial_s - a * alpha - y * d;
            const auto aa = two_third - tensor::stress::double_contraction(d);
            const auto bb = tensor::stress::double_contraction(d, ref);
            const auto cc = tensor::stress::double_contraction(ref);
            const auto sqrt_term = sqrt(bb * bb + aa * cc);

            const auto current_s = tensor::dev(current_stress);
            const vec incre_s = trial_s - current_s;

            const vec base = current_s - a * alpha - y * d;
            const auto incre_incre = tensor::stress::double_contraction(incre_s);
            const auto incre_d = tensor::stress::double_contraction(incre_s, d);

            auto x = .5;
            auto inner_counter = 0u;
            while(true) {
                if(max_iteration == ++inner_counter) {
                    suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
                    return SUANPAN_FAIL;
                }

                const vec middle = base + x * incre_s;
                const auto middle_d = tensor::stress::double_contraction(d, middle);
                const auto tmp_sqrt = std::max(datum::eps, sqrt(middle_d * middle_d + aa * tensor::stress::double_contraction(middle)));
                const auto tmp_numerator = middle_d * incre_d + aa * tensor::stress::double_contraction(incre_s, middle);
                const auto residual_x = tmp_sqrt * incre_d + tmp_numerator;
                const auto jacobian_x = incre_d * residual_x + tmp_sqrt * aa * incre_incre;
                const auto incre_x = tmp_sqrt * residual_x / jacobian_x;

                if(!std::isfinite(incre_x)) {
                    suanpan_error("Infinite number detected.\n");
                    return SUANPAN_FAIL;
                }

                const auto error = fabs(incre_x);
                if(1u == counter) ref_error = error;
                suanpan_debug("Local initial yield ratio iteration error: {:.5E}.\n", error);
                if(error < tolerance * ref_error || ((error < tolerance || fabs(residual_x) < tolerance) && inner_counter > 3u)) {
                    if(x >= 1.) {
                        // elastic unloading
                        z = (bb + sqrt_term) / aa / y;
                        return SUANPAN_SUCCESS;
                    }
                    if(x > 0.) {
                        start_z = (middle_d + tmp_sqrt) / aa / y;
                        suanpan_debug("Initial yield ratio: {:.5E}, corrected yield ratio: {:.5E}.\n", current_z, start_z);
                    }
                    break;
                }

                x -= incre_x;
            }
        }

        const auto trial_ratio = yield_ratio(z);
        const auto avg_rate = u * trial_ratio[0];

        residual(0) = tensor::stress::norm(trial_s - gamma * double_shear * n - a * alpha + (z - 1.) * y * d) - root_two_third * z * y;
        residual(1) = z - start_z - root_two_third * gamma * avg_rate;

        jacobian(0, 0) = tensor::stress::double_contraction(n, pzetapgamma) - double_shear - root_two_third * (b.rb() / bot_alpha * (a + gamma * da - gamma * a * b.r() / bot_alpha) + c.rb() * (1. - z) / bot_d * (y + gamma * dy - gamma * y * c.r() / bot_d) + z * dy);
        jacobian(0, 1) = tensor::stress::double_contraction(n, pzetapz) + root_two_third * y * (gamma * c.rb() / bot_d - 1.);

        jacobian(1, 0) = -root_two_third * avg_rate;
        jacobian(1, 1) = 1. - root_two_third * gamma * u * trial_ratio[1];

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) {
            if(gamma < 0.) {
                suanpan_error("Somehow the plastic multiplier is negative, likely a bug.\n");
                return SUANPAN_FAIL;
            }
            iteration = counter;
            trial_stress -= gamma * double_shear * n;
            trial_stiffness -= double_shear * double_shear * gamma / norm_zeta * unit_dev_tensor - double_shear * double_shear * (gamma / norm_zeta + jacobian(1, 1) / det(jacobian)) * n * n.t();
            return SUANPAN_SUCCESS;
        }

        gamma -= incre(0);
        z -= incre(1);
        if(z > 1.) z = 1. - datum::eps;
        else if(z < 0.) z = 0.;
    }
}

int Subloading::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int Subloading::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int Subloading::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void Subloading::print() {
    suanpan_info("A 3D combined hardening material using subloading surface model.\n");
}
