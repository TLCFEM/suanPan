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

#include "Subloading.h"

#include <Toolbox/tensor.h>

const mat Subloading::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

Subloading::Subloading(const unsigned T, DataSubloading&& D, const double R)
    : DataSubloading{std::move(D)}
    , Material3D(T, R) { access::rw(tolerance) = 1E-13; }

int Subloading::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic, poisson);

    initialize_history(15);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Subloading::get_copy() { return std::make_unique<Subloading>(*this); }

int Subloading::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    const auto norm_incre_strain = tensor::strain::norm(incre_strain);
    if(norm_incre_strain <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    const auto& current_q = current_history(1);
    const auto& current_z = current_history(2);
    const vec current_na(&current_history(3), 6, false, true);
    const vec current_nd(&current_history(9), 6, false, true);
    auto& iteration = trial_history(0);
    auto& q = trial_history(1);
    auto& z = trial_history(2);
    vec na(&trial_history(3), 6, false, true);
    vec nd(&trial_history(9), 6, false, true);

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

        const auto incre_q = root_two_third * gamma;
        q = current_q + incre_q;

        auto [y, dy] = iso_bound(q, true);
        y *= root_two_third;
        dy *= 2. / 3.;

        auto [a, da] = kin_bound(q, true);
        a *= root_two_third;
        da *= 2. / 3.;

        const auto bot_na = 1. + b.r() * incre_q;
        const auto bot_nd = 1. + c.r() * incre_q;

        const vec pzetapz = y / bot_nd * current_nd;
        const vec pzetapg = (b.r() * root_two_third * a / bot_na - da) / bot_na * current_na + (z - 1.) * (dy - c.r() * root_two_third * y / bot_nd) / bot_nd * current_nd;

        const vec zeta = trial_s - a / bot_na * current_na + (z - 1.) * pzetapz;
        const auto norm_zeta = tensor::stress::norm(zeta);
        const vec n = zeta / norm_zeta;

        na = (incre_q * b.rb() * n + current_na) / bot_na;
        nd = (incre_q * c.rb() * n + current_nd) / bot_nd;

        if(1u == counter) {
            const vec ref = trial_s - a * na - y * nd;
            const vec incre_s = trial_s - tensor::dev(current_stress);
            const vec base = ref - incre_s;

            const auto aa = 1. - tensor::stress::double_contraction(nd);
            const auto bb = tensor::stress::double_contraction(nd, ref);
            const auto cc = tensor::stress::double_contraction(ref);

            const auto incre_incre = tensor::stress::double_contraction(incre_s);
            const auto incre_d = tensor::stress::double_contraction(incre_s, nd);

            auto x = .5;
            auto inner_counter = 0u;
            while(true) {
                if(max_iteration == ++inner_counter) {
                    suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
                    return SUANPAN_FAIL;
                }

                const vec middle = base + x * incre_s;
                const auto middle_d = tensor::stress::double_contraction(nd, middle);
                const auto tmp_sqrt = std::max(datum::eps, std::sqrt(middle_d * middle_d + aa * tensor::stress::double_contraction(middle)));
                const auto residual_x = tmp_sqrt * incre_d + middle_d * incre_d + aa * tensor::stress::double_contraction(incre_s, middle);
                const auto jacobian_x = incre_d * residual_x + tmp_sqrt * aa * incre_incre;
                const auto incre_x = tmp_sqrt * residual_x / jacobian_x;

                if(!std::isfinite(incre_x)) {
                    suanpan_error("Infinite number detected.\n");
                    return SUANPAN_FAIL;
                }

                const auto error = std::fabs(incre_x);
                if(1u == inner_counter) ref_error = error;
                suanpan_debug("Local initial yield ratio iteration error: {:.5E}.\n", error);
                if(error < tolerance * ref_error || ((error < tolerance || std::fabs(residual_x) < tolerance) && inner_counter > 3u)) {
                    if(x >= 1.) {
                        // elastic unloading
                        z = (bb + std::sqrt(bb * bb + aa * cc)) / aa / y;
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

        residual(0) = tensor::stress::norm(trial_s - gamma * double_shear * n - a * na + (z - 1.) * y * nd) - z * y;
        residual(1) = z - start_z - incre_q * avg_rate;

        jacobian(0, 0) = tensor::stress::double_contraction(n, pzetapg) - double_shear - root_two_third * (b.rb() / bot_na * (a + gamma * da - incre_q * a / bot_na * b.r()) + (1. - z) * c.rb() / bot_nd * (y + gamma * dy - incre_q * y / bot_nd * c.r())) - z * dy;
        jacobian(0, 1) = tensor::stress::double_contraction(n, pzetapz) + y * incre_q * c.rb() / bot_nd - y;

        jacobian(1, 0) = -root_two_third * avg_rate;
        jacobian(1, 1) = 1. - incre_q * u * trial_ratio[1];

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

        const auto error = suanpan::inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || suanpan::inf_norm(residual) < tolerance) && counter > 5u)) {
            if(gamma < 0.) {
                suanpan_error("Somehow the plastic multiplier is negative, likely a bug.\n");
                return SUANPAN_FAIL;
            }
            iteration = counter;
            trial_stress -= gamma * double_shear * n;
            trial_stiffness -= double_shear * double_shear * gamma / norm_zeta * unit_dev_tensor - double_shear * double_shear * (gamma / norm_zeta + jacobian(1, 1) / det(jacobian)) * n * n.t();
            return SUANPAN_SUCCESS;
        }

        gamma = suanpan::clamp(gamma - incre(0), 0., 1.1 * norm_incre_strain);
        z = suanpan::clamp_unit(z - incre(1));
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
    suanpan_info("A 3D combined hardening material using subloading surface model. doi:10.1007/s00707-025-04339-0\n");
}
