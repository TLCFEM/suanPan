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

#include "Balloon1D.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

const double DataBalloon1D::Saturation::root_one_half = sqrt(1.5);

const double Balloon1D::rate_bound = -log(z_bound);

pod2 Balloon1D::yield_ratio(const double z) {
    if(z < z_bound) return {rate_bound, 0.};

    return {-log(z), -1. / z};
}

Balloon1D::Balloon1D(const unsigned T, DataBalloon1D&& D, const double R)
    : DataBalloon1D{std::move(D)}
    , Material1D(T, R) {}

int Balloon1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic;

    initialize_history(6u + static_cast<unsigned>(b.size() + c.size()));

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Balloon1D::get_copy() { return std::make_unique<Balloon1D>(*this); }

int Balloon1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(std::fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    trial_zr = current_zr;
    const auto& current_q = current_history(2);
    const auto& current_qm = current_history(3);
    const auto& current_qr = current_history(4);
    const auto& current_z = current_history(5);
    auto& iteration = trial_history(0);
    auto& last_loading = trial_history(1);
    auto& q = trial_history(2);
    auto& qm = trial_history(3);
    auto& qr = trial_history(4);
    auto& z = trial_history(5);

    const vec current_alpha(&current_history(6), b.size(), false, true);
    const vec current_d(&current_history(6 + b.size()), c.size(), false, true);
    vec alpha(&trial_history(6), b.size(), false, true);
    vec d(&trial_history(6 + b.size()), c.size(), false, true);

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

        auto split = 1., dsplit = 0.;
        if(const auto max_zr = trial_zr.max(); z < max_zr) {
            if(k >= 1.) split = 0.;
            else {
                const auto x = z / max_zr;
                const auto denom = 1. + k - 2. * k * x;
                split = (1. - k) / denom * x;
                dsplit = (1. - k) / denom * (1. + k) / denom;
            }
        }

        const auto [ym, dym] = iso_m(qm = current_qm + split * gamma);
        const auto [yr, dyr] = iso_r(qr = current_qr + (1. - split) * gamma);

        const auto y = ym + yr;
        const auto pypg = (dym - dyr) * split + dyr;
        const auto pypz = (dym - dyr) * gamma * dsplit;

        const auto [a, da] = kin(q = current_q + gamma);

        vec bottom_alpha(b.size(), fill::none), bottom_d(c.size(), fill::none);
        for(auto I = 0llu; I < b.size(); ++I) bottom_alpha(I) = 1. + b[I].r() * gamma;
        for(auto I = 0llu; I < c.size(); ++I) bottom_d(I) = 1. + c[I].r() * gamma;

        const auto n = trial_stress(0) - a * accu(current_alpha / bottom_alpha) + (z - 1.) * y * accu(current_d / bottom_d) > 0. ? 1. : -1.;

        for(auto I = 0llu; I < b.size(); ++I) alpha(I) = (b[I].rb() * gamma * n + current_alpha(I)) / bottom_alpha(I);
        for(auto I = 0llu; I < c.size(); ++I) d(I) = (c[I].rb() * gamma * n + current_d(I)) / bottom_d(I);

        const auto sum_alpha = accu(alpha), sum_d = accu(d);

        if(1u == counter) {
            const auto s = (y * sum_d + a * sum_alpha - current_stress(0)) / (trial_stress(0) - current_stress(0));
            if(std::signbit(last_loading) == std::signbit(s)) trial_zr.enqueue(z);
            if(s >= 1.) {
                last_loading = -1.;
                // elastic unloading
                z = ((trial_stress(0) - a * sum_alpha) / y - sum_d) / (n - sum_d);
                return SUANPAN_SUCCESS;
            }
            last_loading = 1.;
            if(s > 0.) start_z = 0.;
        }

        auto dalpha = 0., dd = 0.;
        for(auto I = 0llu; I < b.size(); ++I) dalpha += b[I].r() * (b[I].b() * n - alpha[I]) / bottom_alpha[I];
        for(auto I = 0llu; I < c.size(); ++I) dd += c[I].r() * (c[I].b() * n - d[I]) / bottom_d[I];

        const auto trial_ratio = yield_ratio(z);
        const auto avg_rate = u * trial_ratio[0];

        residual(0) = std::fabs(trial_stress(0) - elastic * gamma * n - a * sum_alpha + (z - 1.) * y * sum_d) - z * y;
        residual(1) = z - start_z - gamma * avg_rate;

        jacobian(0, 0) = n * ((z - 1.) * (y * dd + sum_d * pypg) - (a * dalpha + sum_alpha * da)) - elastic - z * pypg;
        jacobian(0, 1) = n * sum_d * (y + z * pypz - pypz) - y - z * pypz;

        jacobian(1, 0) = -avg_rate;
        jacobian(1, 1) = 1. - u * gamma * trial_ratio[1];

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) {
            iteration = counter;
            trial_stress -= elastic * gamma * n;
            trial_stiffness += elastic / det(jacobian) * elastic * jacobian(1, 1);
            return SUANPAN_SUCCESS;
        }

        gamma -= incre(0);
        if(gamma < 0.) gamma = 0.;

        z = suanpan::clamp_unit(z - incre(1));
    }
}

int Balloon1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_zr.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int Balloon1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_zr = trial_zr;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int Balloon1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_zr = current_zr;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void Balloon1D::print() {
    suanpan_info("A uniaxial combined hardening material using subloading surface model. doi:10.1007/s00707-025-04339-0\n");
    Material1D::print();
}
