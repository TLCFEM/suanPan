/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "SubloadingMetal.h"
#include <Toolbox/tensor.h>

const double SubloadingMetal::root_two_third = sqrt(two_third);
const double SubloadingMetal::rate_bound = -log(z_bound);
const mat SubloadingMetal::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

vec2 SubloadingMetal::yield_ratio(const double z) {
    if(z < z_bound) return {rate_bound, 0.};

    return {-log(z), -1. / z};
}

SubloadingMetal::SubloadingMetal(const unsigned T, DataSubloadingMetal&& D, const double R)
    : DataSubloadingMetal{std::move(D)}
    , Material3D(T, R) {}

int SubloadingMetal::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic, poissons_ratio);

    initialize_history(14);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> SubloadingMetal::get_copy() { return make_unique<SubloadingMetal>(*this); }

int SubloadingMetal::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    const auto& current_q = current_history(0);
    const auto& current_z = current_history(1);
    const vec current_alpha(&current_history(2), 6, false, true);
    const vec current_d(&current_history(8), 6, false, true);
    auto& q = trial_history(0);
    auto& z = trial_history(1);
    vec alpha(&trial_history(2), 6, false, true);
    vec d(&trial_history(8), 6, false, true);

    const vec trial_s = tensor::dev(trial_stress);

    auto gamma = 0., ref_error = 0.;

    vec2 residual, incre;
    mat22 jacobian;

    const auto current_ratio = yield_ratio(current_z);

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

        const auto bot_alpha = 1. + be * gamma;
        const auto bot_d = 1. + ce * gamma;

        const vec pzetapz = y / bot_d * current_d;
        const vec pzetapgamma = (be * a / bot_alpha - da) / bot_alpha * current_alpha + (z - 1.) * (dy - ce * y / bot_d) / bot_d * current_d;

        const vec zeta = trial_s - a / bot_alpha * current_alpha + (z - 1.) * pzetapz;
        const auto norm_zeta = tensor::stress::norm(zeta);
        const vec n = zeta / norm_zeta;

        alpha = (root_two_third * gamma * be * n + current_alpha) / bot_alpha;

        d = (root_two_third * gamma * ce * ze * n + current_d) / bot_d;

        const vec eta = trial_s - gamma * double_shear * n - a * alpha + (z - 1.) * y * d;

        residual(0) = tensor::stress::norm(eta) - root_two_third * z * y;

        if(1u == counter && residual(0) < 0.) {
            const vec ref = trial_s - a * alpha - y * d;
            const auto aa = (tensor::stress::double_contraction(d) - two_third) * y * y;
            const auto bb = -y * tensor::stress::double_contraction(d, ref);
            const auto cc = tensor::stress::double_contraction(ref);
            const auto sqrt_term = sqrt(bb * bb - aa * cc);

            z = (bb - sqrt_term) / aa;
            if(z < 0. || z > 1.) z = (bb + sqrt_term) / aa;

            return SUANPAN_SUCCESS;
        }

        const auto trial_ratio = yield_ratio(z);
        const auto avg_rate = u * .5 * (current_ratio(0) + trial_ratio(0));

        residual(1) = z - current_z - root_two_third * gamma * avg_rate;

        jacobian(0, 0) = tensor::stress::double_contraction(n, pzetapgamma) - double_shear - root_two_third * (be / bot_alpha * (a + gamma * da - gamma * a * be / bot_alpha) + ce * ze * (1. - z) / bot_d * (y + gamma * dy - gamma * y * ce / bot_d) + z * dy);
        jacobian(0, 1) = tensor::stress::double_contraction(n, pzetapz) + root_two_third * y * (gamma * ce * ze / bot_d - 1.);

        jacobian(1, 0) = -root_two_third * avg_rate;
        jacobian(1, 1) = 1. - root_two_third * gamma * u * .5 * trial_ratio(1);

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) {
            trial_stress -= gamma * double_shear * n;

            mat right(2, 6, fill::zeros);
            right.row(0) = -n.t();

            const mat left = solve(jacobian, right);

            trial_stiffness -= double_shear * double_shear * n * left.row(0) + double_shear * double_shear * gamma / norm_zeta * (unit_dev_tensor - n * n.t());

            return SUANPAN_SUCCESS;
        }

        gamma -= incre(0);
        z -= incre(1);
        if(z > 1.) z = 1. - datum::eps;
    }
}

int SubloadingMetal::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int SubloadingMetal::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int SubloadingMetal::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void SubloadingMetal::print() {
    suanpan_info("A 3D combined hardening material using subloading surface model.\n");
}
