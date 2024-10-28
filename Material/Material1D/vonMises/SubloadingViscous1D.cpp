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

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    const auto& current_q = current_history(1);
    const auto& current_z = current_history(2);
    const auto& current_zv = current_history(3);
    auto& iteration = trial_history(0);
    auto& q = trial_history(1);
    auto& z = trial_history(2);
    auto& zv = trial_history(3);

    const vec current_alpha(&current_history(4), b.size(), false, true);
    const vec current_d(&current_history(4 + b.size()), c.size(), false, true);
    vec alpha(&trial_history(4), b.size(), false, true);
    vec d(&trial_history(4 + b.size()), c.size(), false, true);

    const auto incre_t = *incre_time > 0. ? *incre_time : 1.;

    iteration = 0.;
    auto gamma = 0., ref_error = 0.;
    auto start_z = current_z;

    vec3 residual, incre;
    mat33 jacobian;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        q = current_q + gamma;

        const auto exp_iso = saturation_iso * exp(-m_iso * q);
        auto y = initial_iso + saturation_iso + k_iso * q - exp_iso;
        auto dy = k_iso + m_iso * exp_iso;
        if(y < 0.) y = dy = 0.;

        const auto exp_kin = saturation_kin * exp(-m_kin * q);
        auto a = initial_kin + saturation_kin + k_kin * q - exp_kin;
        auto da = k_kin + m_kin * exp_kin;
        if(a < 0.) a = da = 0.;

        vec bottom_alpha(b.size(), fill::none), bottom_d(c.size(), fill::none);
        for(auto I = 0llu; I < b.size(); ++I) bottom_alpha(I) = 1. + b[I].r() * gamma;
        for(auto I = 0llu; I < c.size(); ++I) bottom_d(I) = 1. + c[I].r() * gamma;

        const auto n = trial_stress(0) - a * sum(current_alpha / bottom_alpha) + (zv - 1.) * y * sum(current_d / bottom_d) > 0. ? 1. : -1.;

        for(auto I = 0llu; I < b.size(); ++I) alpha(I) = (b[I].rb() * gamma * n + current_alpha(I)) / bottom_alpha(I);
        for(auto I = 0llu; I < c.size(); ++I) d(I) = (c[I].rb() * gamma * n + current_d(I)) / bottom_d(I);

        const auto sum_alpha = sum(alpha), sum_d = sum(d);

        if(1u == counter) {
            const auto s = (y * sum_d + a * sum_alpha - current_stress(0)) / (trial_stress(0) - current_stress(0));
            if(s >= 1.) {
                // elastic unloading
                zv = ((trial_stress(0) - a * sum_alpha) / y - sum_d) / (n - sum_d);
                z = zv / current_zv * current_z;
                return SUANPAN_SUCCESS;
            }
            if(s > 0.) start_z = 0.;
        }

        auto dalpha = 0., dd = 0.;
        for(auto I = 0llu; I < b.size(); ++I) dalpha += b[I].r() * (b[I].b() * n - alpha[I]) / bottom_alpha[I];
        for(auto I = 0llu; I < c.size(); ++I) dd += c[I].r() * (c[I].b() * n - d[I]) / bottom_d[I];

        const auto trial_ratio = yield_ratio(z);
        const auto avg_rate = u * trial_ratio(0);

        residual(0) = fabs(trial_stress(0) - elastic * gamma * n - a * sum_alpha + (zv - 1.) * y * sum_d) - zv * y;
        residual(1) = z - start_z - gamma * avg_rate;
        residual(2) = (zv - cv) * mu * gamma;

        jacobian(0, 0) = n * ((zv - 1.) * (y * dd + sum_d * dy) - (a * dalpha + sum_alpha * da)) - elastic - zv * dy;
        jacobian(0, 1) = n * y * sum_d - y;
        jacobian(0, 2) = 0.;

        jacobian(1, 0) = -avg_rate;
        jacobian(1, 1) = 0.;
        jacobian(1, 2) = 1. - u * gamma * trial_ratio(1);

        jacobian(2, 0) = (zv - cv) * mu;
        jacobian(2, 1) = mu * gamma;
        jacobian(2, 2) = 0.;

        if(zv > z) {
            const auto diff_z = zv - z;
            const auto power_term = incre_t * pow(diff_z, nv - 1);

            residual(2) += diff_z * power_term;

            jacobian(2, 1) += nv * power_term;
            jacobian(2, 2) -= nv * power_term;
        }

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
            trial_stress -= elastic * gamma * n;
            trial_stiffness += elastic / det(jacobian) * elastic * det(jacobian.submat(1, 1, 2, 2));
            return SUANPAN_SUCCESS;
        }

        gamma -= incre(0);
        zv -= incre(1);
        z -= incre(2);
        if(z > 1.) z = 1. - datum::eps;
        else if(z < 0.) z = 0.;
    }
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
    suanpan_info("A uniaxial combined hardening material using subloading surface model.\n");
    Material1D::print();
}
