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

#include "SubloadingViscous1D.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

const double DataSubloadingViscous1D::Saturation::root_one_half = sqrt(1.5);

const double SubloadingViscous1D::rate_bound = -log(z_bound);

vec2 SubloadingViscous1D::yield_ratio(const double z) {
    if(z < z_bound) return {rate_bound, 0.};

    return {-log(z), -1. / z};
}

std::tuple<double, double> SubloadingViscous1D::isotropic_bound(const double q) const {
    const auto exp_iso = saturation_iso * exp(-m_iso * q);
    auto y = initial_iso + saturation_iso + k_iso * q - exp_iso;
    auto dy = k_iso + m_iso * exp_iso;
    if(y < 0.) y = dy = 0.;
    return {y, dy};
}

std::tuple<double, double> SubloadingViscous1D::kinematic_bound(const double q) const {
    const auto exp_kin = saturation_kin * exp(-m_kin * q);
    auto a = initial_kin + saturation_kin + k_kin * q - exp_kin;
    auto da = k_kin + m_kin * exp_kin;
    if(a < 0.) a = da = 0.;
    return {a, da};
}

int SubloadingViscous1D::partial_loading(double& fragment, vec& start_history, const double start_stress, const double diff_stress) {
    const auto& current_q = start_history(1);
    const auto& current_z = start_history(2);
    const vec current_alpha(&start_history(4), b.size(), false, true);
    const vec current_d(&start_history(4 + b.size()), c.size(), false, true);

    trial_history = start_history;
    auto& iteration = trial_history(0);
    auto& q = trial_history(1);
    auto& z = trial_history(2);
    auto& zv = trial_history(3);
    vec alpha(&trial_history(4), b.size(), false, true);
    vec d(&trial_history(4 + b.size()), c.size(), false, true);

    const auto norm_mu = mu / (incre_time && *incre_time > 0. ? *incre_time : 1.);

    iteration = 0.;
    auto gamma = 0., ref_error = 0.;

    vec4 residual, incre;
    mat44 jacobian;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        if(fragment > 1.) fragment = 1.;

        residual.zeros();
        jacobian.zeros();

        const auto partial_stress = start_stress + fragment * diff_stress;

        q = current_q + gamma;

        const auto [y, dy] = isotropic_bound(q);
        const auto [a, da] = kinematic_bound(q);

        vec bottom_alpha(b.size(), fill::none), bottom_d(c.size(), fill::none);
        for(auto I = 0llu; I < b.size(); ++I) bottom_alpha(I) = 1. + b[I].r() * gamma;
        for(auto I = 0llu; I < c.size(); ++I) bottom_d(I) = 1. + c[I].r() * gamma;

        const auto n = partial_stress - a * accu(current_alpha / bottom_alpha) + (z - 1.) * y * accu(current_d / bottom_d) > 0. ? 1. : -1.;

        for(auto I = 0llu; I < b.size(); ++I) alpha(I) = (b[I].rb() * gamma * n + current_alpha(I)) / bottom_alpha(I);
        for(auto I = 0llu; I < c.size(); ++I) d(I) = (c[I].rb() * gamma * n + current_d(I)) / bottom_d(I);

        const auto sum_alpha = accu(alpha), sum_d = accu(d);

        auto dalpha = 0., dd = 0.;
        for(auto I = 0llu; I < b.size(); ++I) dalpha += b[I].r() * (b[I].b() * n - alpha[I]) / bottom_alpha[I];
        for(auto I = 0llu; I < c.size(); ++I) dd += c[I].r() * (c[I].b() * n - d[I]) / bottom_d[I];

        const auto trial_ratio = yield_ratio(z);
        const auto avg_rate = u * trial_ratio(0);
        const auto fraction_term = (cv * z - zv) * norm_mu * gamma + 1.;
        const auto power_term = pow(fraction_term, nv - 1.);

        residual(0) = fabs(partial_stress - elastic * gamma * n - a * sum_alpha + (z - 1.) * y * sum_d) - zv * y;
        residual(1) = z - current_z - gamma * avg_rate;
        residual(2) = zv - fraction_term * power_term * z;

        jacobian(0, 0) = n * ((z - 1.) * (y * dd + sum_d * dy) - (a * dalpha + sum_alpha * da)) - elastic - zv * dy;
        jacobian(0, 1) = -y;
        jacobian(0, 2) = n * y * sum_d;

        jacobian(1, 0) = -avg_rate;
        jacobian(1, 2) = 1. - u * gamma * trial_ratio(1);

        jacobian(2, 0) = -z * nv * power_term * (cv * z - zv) * norm_mu;
        jacobian(2, 1) = 1. + z * nv * power_term * norm_mu * gamma;
        jacobian(2, 2) = -power_term * (fraction_term + z * nv * cv * norm_mu * gamma);

        if(fragment >= 1.) jacobian(3, 3) = 1.;
        else {
            residual(3) = zv - z;

            jacobian(3, 1) = 1.;
            jacobian(3, 2) = -1.;

            jacobian(0, 3) = n * diff_stress;
        }

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) {
            iteration = counter;
            trial_stress = partial_stress - elastic * gamma * n;
            if(fragment >= 1.) trial_stiffness += elastic / det(jacobian.submat(0, 0, 2, 2)) * elastic * det(jacobian.submat(1, 1, 2, 2));
            return SUANPAN_SUCCESS;
        }

        gamma -= incre(0);
        zv -= incre(1);
        z -= incre(2);
        fragment -= incre(3);
        if(gamma < 0.) gamma = 0.;
        if(z < 0.) z = 0.;
        else if(z > 1.) z = 1.;
        if(zv < z) zv = z;
        else if(zv > cv * z) zv = cv * z;
    }
}

SubloadingViscous1D::SubloadingViscous1D(const unsigned T, DataSubloadingViscous1D&& D, const double R)
    : DataSubloadingViscous1D{std::move(D)}
    , Material1D(T, R) {}

int SubloadingViscous1D::initialize(const shared_ptr<DomainBase>& D) {
    if(nullptr != D) incre_time = &D->get_factory()->modify_incre_time();

    trial_stiffness = current_stiffness = initial_stiffness = elastic;

    initialize_history(4u + static_cast<unsigned>(b.size() + c.size()));

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> SubloadingViscous1D::get_copy() { return make_unique<SubloadingViscous1D>(*this); }

int SubloadingViscous1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;
    incre_stress = (trial_stiffness = initial_stiffness) * incre_strain;
    trial_stress = current_stress + incre_stress;

    const auto& current_q = current_history(1);
    const auto& current_z = current_history(2);
    const auto current_sum_alpha = sum(vec(&current_history(4), b.size(), false, true));
    const auto current_sum_d = sum(vec(&current_history(4 + b.size()), c.size(), false, true));

    const auto [current_y, current_dy] = isotropic_bound(current_q);
    const auto [current_a, current_da] = kinematic_bound(current_q);

    auto fragment = 0.;

    // pure plastic loading
    if(const auto current_eta = current_stress(0) - current_a * current_sum_alpha + (current_z - 1.) * current_y * current_sum_d; current_eta * incre_stress(0) >= 0.) return partial_loading(fragment = 1., current_history, trial_stress(0), 0.);

    if(SUANPAN_SUCCESS != partial_loading(fragment, current_history, current_stress(0), incre_stress(0))) return SUANPAN_FAIL;

    // pure plastic unloading
    if(fragment >= 1.) return SUANPAN_SUCCESS;

    auto inter_history = trial_history;
    const auto& inter_q = inter_history(1);
    auto& inter_z = inter_history(2);
    auto& inter_zv = inter_history(3);
    const auto inter_sum_alpha = sum(vec(&inter_history(4), b.size(), false, true));
    const auto inter_sum_d = sum(vec(&inter_history(4 + b.size()), c.size(), false, true));

    const auto [inter_y, inter_dy] = isotropic_bound(inter_q);
    const auto [inter_a, inter_da] = kinematic_bound(inter_q);

    const auto remaining_stress = (1. - fragment) * incre_stress(0);
    const auto inter_stress = trial_stress(0);
    trial_stress += remaining_stress;

    if(const auto inter_centre = inter_a * inter_sum_alpha + inter_y * inter_sum_d; remaining_stress > 0. ? inter_centre >= trial_stress(0) : inter_centre <= trial_stress(0)) {
        // pure elastic unloading
        const auto inter_n = inter_stress - inter_a * inter_sum_alpha + (inter_z - 1.) * inter_y * inter_sum_d > 0. ? 1. : -1.;
        trial_history(3) = trial_history(2) = ((trial_stress(0) - inter_a * inter_sum_alpha) / inter_y - inter_sum_d) / (inter_n - inter_sum_d);
        return SUANPAN_SUCCESS;
    }

    // plastic loading after unloading
    inter_zv = inter_z = 0.;

    return partial_loading(fragment = 1., inter_history, trial_stress(0), 0.);

    // auto& iteration = trial_history(0);
    // auto& q = trial_history(1);
    // auto& z = trial_history(2);
    // auto& zv = trial_history(3);
    // vec alpha(&trial_history(4), b.size(), false, true);
    // vec d(&trial_history(4 + b.size()), c.size(), false, true);
}

int SubloadingViscous1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int SubloadingViscous1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int SubloadingViscous1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void SubloadingViscous1D::print() {
    suanpan_info("A uniaxial combined hardening material using subloading surface model with optional viscosity.\n");
    Material1D::print();
}
