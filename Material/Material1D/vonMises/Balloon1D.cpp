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

Balloon1D::Balloon1D(const unsigned T, DataBalloon1D&& D, const double R)
    : DataBalloon1D{std::move(D)}
    , Material1D(T, R) {}

int Balloon1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic;

    initialize_history(6u + static_cast<unsigned>(ba.size() + bb.size() + bd.size()));

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
    const auto& current_qc = current_history(4);
    const auto& current_z = current_history(5);
    auto& iteration = trial_history(0);
    auto& last_loading = trial_history(1);
    auto& q = trial_history(2);
    auto& qm = trial_history(3);
    auto& qc = trial_history(4);
    auto& z = trial_history(5);

    const vec current_alpha(&current_history(6), ba.size(), false, true);
    const vec current_beta(&current_history(6 + ba.size()), bb.size(), false, true);
    const vec current_d(&current_history(6 + ba.size() + bb.size()), bd.size(), false, true);
    vec alpha(&trial_history(6), ba.size(), false, true);
    vec beta(&trial_history(6 + ba.size()), bb.size(), false, true);
    vec d(&trial_history(6 + ba.size() + bb.size()), bd.size(), false, true);

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

        auto km = 1.;
        if(const auto ref_zr = zr_size > 0 ? trial_zr.max() : trial_zr.min(); z < ref_zr) km = kb / (1. - kr * std::log(1. - ref_zr));
        const auto kc = 1. - km;

        qm = current_qm + km * gamma;
        qc = current_qc + kc * gamma;
        q = current_q + gamma;

        const auto [hf, dhf] = bound_hf(qm, true);
        const auto phfpg = dhf * km;

        const auto [ha, dha] = bound_ha(q, true);
        const auto phapg = dha;

        const auto [hb, dhb] = bound_hb(qm, false);
        const auto phbpg = dhb * km;

        const auto [hd, dhd] = bound_hd(q, true);
        const auto phdpg = dhd;

        vec top_alpha(ba.size(), fill::none), top_beta(bb.size(), fill::none), top_d(bd.size(), fill::none);
        vec bot_alpha(ba.size(), fill::none), bot_beta(bb.size(), fill::none), bot_d(bd.size(), fill::none);
        for(auto I = 0llu; I < ba.size(); ++I) {
            top_alpha(I) = ba[I].rb() * gamma;
            bot_alpha(I) = 1. + ba[I].r() * gamma;
        }
        for(auto I = 0llu; I < bb.size(); ++I) {
            top_beta(I) = bb[I].rb() * gamma * kc;
            bot_beta(I) = 1. + bb[I].r() * gamma * kc;
        }
        for(auto I = 0llu; I < bd.size(); ++I) {
            top_d(I) = bd[I].rb() * gamma;
            bot_d(I) = 1. + bd[I].r() * gamma;
        }

        const auto zeta = trial_stress(0) - ha * accu(current_alpha / bot_alpha) - hb * accu(current_beta / bot_beta) + (z - 1.) * hd * accu(current_d / bot_d);
        const auto norm_zeta = std::fabs(zeta);
        const auto norm_diff = elastic * gamma + ha * accu(top_alpha / bot_alpha) + hb * accu(top_beta / bot_beta) + (1. - z) * hd * accu(top_d / bot_d);

        auto n = 0.;
        if(norm_zeta >= norm_diff) n = zeta > 0. ? 1. : -1.;
        else if(-norm_zeta >= norm_diff) n = zeta > 0. ? -1. : 1.;
        else {
            suanpan_error("Sign mismatch, likely an unknown bug.\n");
            return SUANPAN_FAIL;
        }

        const auto sum_alpha = accu(alpha = (top_alpha * n + current_alpha) / bot_alpha);
        const auto sum_beta = accu(beta = (top_beta * n + current_beta) / bot_beta);
        const auto sum_d = accu(d = (top_d * n + current_d) / bot_d);

        if(1u == counter) {
            const auto h_sum = ha * sum_alpha + hb * sum_beta + hd * sum_d;
            const auto s = (h_sum - current_stress(0)) / (trial_stress(0) - current_stress(0));
            if(std::signbit(last_loading) == std::signbit(s)) trial_zr.enqueue(z);
            if(s >= 1.) {
                last_loading = -1.;
                // elastic unloading
                z = (trial_stress(0) - h_sum) / (hf * n - hd * sum_d);
                return SUANPAN_SUCCESS;
            }
            last_loading = 1.;
            if(s > 0.) start_z = 0.;
        }

        auto dalpha = 0., dbeta = 0., dd = 0.;
        for(auto I = 0llu; I < ba.size(); ++I) dalpha += ba[I].r() * (ba[I].b() * n - alpha[I]) / bot_alpha[I];
        for(auto I = 0llu; I < bb.size(); ++I) dbeta += bb[I].r() * kc * (bb[I].b() * n - beta[I]) / bot_beta[I];
        for(auto I = 0llu; I < bd.size(); ++I) dd += bd[I].r() * (bd[I].b() * n - d[I]) / bot_d[I];

        const auto trial_ratio = yield_ratio(z);
        const auto avg_rate = u * trial_ratio[0];
        const auto diff_z = z - start_z;

        residual(0) = std::fabs(trial_stress(0) - elastic * gamma * n - ha * sum_alpha - hb * sum_beta + (z - 1.) * hd * sum_d) - z * hf;
        residual(1) = hf * diff_z - gamma * avg_rate;

        jacobian(0, 0) = n * ((z - 1.) * (hd * dd + sum_d * phdpg) - ha * dalpha - sum_alpha * phapg - hb * dbeta - sum_beta * phbpg) - elastic - z * phfpg;
        jacobian(0, 1) = n * hd * sum_d - hf;

        jacobian(1, 0) = phfpg * diff_z - avg_rate;
        jacobian(1, 1) = hf - u * gamma * trial_ratio[1];

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

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
    suanpan_info("The Balloon uniaxial model.\n");
    Material1D::print();
}
