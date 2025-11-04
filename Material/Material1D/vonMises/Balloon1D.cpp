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

/**
 * @brief Perform the initial check to determine whether the material is plastic loading or elastic unloading.
 * If it is plastic loading, correct the $z$ value when necessary.
 * If it is elastic unloading, update the $z$ value accordingly.
 * The loading flag is also updated.
 * The history buffer for $z$ is also updated.
 * @param start_z The $z$ value at the start of the step
 * @return The revised $z$ value at the start of the integration
 */
double Balloon1D::initial_check(double start_z) {
    const auto& qm = current_history(3);

    const vec fc(&current_history(5), bfc.size(), false, true);
    const vec ac(&current_history(5 + bfc.size()), bac.size(), false, true);
    const vec alpha(&current_history(5 + bfc.size() + bac.size()), ba.size(), false, true);
    const vec d(&current_history(5 + bfc.size() + bac.size() + ba.size()), bd.size(), false, true);

    [[maybe_unused]] const auto [fm, dfm] = bound_fm(qm, true);
    [[maybe_unused]] const auto [am, dam] = bound_am(qm, true);

    const auto hf = std::max(0., fm + accu(fc));
    const auto ha = std::max(0., am + accu(ac));

    const auto sum_alpha = ha * accu(alpha);
    const auto sum_d = hf * accu(d);
    const auto sum_all = sum_alpha + sum_d;

    auto& last_loading = trial_history(1);
    auto& z = trial_history(4);

    // assuming z is zero find the load factor
    const auto s = (sum_all - current_stress(0)) / incre_strain(0);
    if(std::signbit(last_loading) == std::signbit(s)) trial_zr.enqueue(z);
    if(s >= elastic) {
        last_loading = -1.;
        // elastic unloading
        const auto n = current_stress(0) + z * sum_d > sum_all ? 1. : -1.;
        z = (trial_stress(0) - sum_all) / (hf * n - sum_d);
    }
    else {
        last_loading = 1.;
        if(s > 0.) start_z = 0.;
    }

    return start_z;
}

auto Balloon1D::compute_isotropic_bound(const double gamma, const double km, const double dkm) {
    const auto& qm = trial_history(3);

    const vec current_hfc(&current_history(5), bfc.size(), false, true);
    vec hfc(&trial_history(5), bfc.size(), false, true);

    const auto kc = 1. - km, dkc = -dkm;

    const auto [fm, dfm] = bound_fm(qm, true);
    const auto [fc, dfc] = bound_fc(qm, false);
    const auto pfcpg = dfc * km, pfcpz = dfc * gamma * dkm;

    auto phfpg = dfm * km, phfpz = dfm * gamma * dkm;
    for(auto I = 0llu; I < bfc.size(); ++I) {
        const auto bot_fc = 1. + bfc[I].b() * gamma * kc;
        hfc(I) = (bfc[I].a() * gamma * kc * fc + current_hfc(I)) / bot_fc;
        phfpg += (bfc[I].a() * (fc + gamma * pfcpg) - hfc(I) * bfc[I].b()) * kc / bot_fc;
        phfpz += (bfc[I].a() * (dkc * fc + kc * pfcpz) - hfc(I) * bfc[I].b() * dkc) * gamma / bot_fc;
    }

    if(const auto hf = fm + accu(hfc); hf > 0.) return std::make_tuple(hf, phfpg, phfpz);

    return std::make_tuple(0., 0., 0.);
}

auto Balloon1D::compute_kinematic_bound(const double gamma, const double km, const double dkm) {
    const auto& qm = trial_history(3);

    const vec current_hac(&current_history(5 + bfc.size()), bac.size(), false, true);
    vec hac(&trial_history(5 + bfc.size()), bac.size(), false, true);

    const auto kc = 1. - km, dkc = -dkm;

    const auto [am, dam] = bound_am(qm, true);
    const auto [ac, dac] = bound_ac(qm, false);
    const auto pacpg = dac * km, pacpz = dac * gamma * dkm;

    auto phapg = dam * km, phapz = dam * gamma * dkm;
    for(auto I = 0llu; I < bac.size(); ++I) {
        const auto bot_ac = 1. + bac[I].b() * gamma * kc;
        hac(I) = (bac[I].a() * gamma * kc * ac + current_hac(I)) / bot_ac;
        phapg += (bac[I].a() * (ac + gamma * pacpg) - hac(I) * bac[I].b()) * kc / bot_ac;
        phapz += (bac[I].a() * (dkc * ac + kc * pacpz) - hac(I) * bac[I].b() * dkc) * gamma / bot_ac;
    }

    if(const auto ha = am + accu(hac); ha > 0.) return std::make_tuple(ha, phapg, phapz);

    return std::make_tuple(0., 0., 0.);
}

Balloon1D::Balloon1D(const unsigned T, DataBalloon1D&& D, const double R)
    : DataBalloon1D{std::move(D)}
    , Material1D(T, R) {}

int Balloon1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic;

    initialize_history(5u + static_cast<unsigned>(bfc.size() + bac.size() + ba.size() + bd.size()));

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Balloon1D::get_copy() { return std::make_unique<Balloon1D>(*this); }

int Balloon1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    const auto abs_incre_strain = std::fabs(incre_strain(0));
    if(abs_incre_strain <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    trial_zr = current_zr;
    const auto& current_q = current_history(2);
    const auto& current_qm = current_history(3);
    const auto& current_z = current_history(4);
    auto& iteration = trial_history(0);
    auto& last_loading = trial_history(1);
    auto& q = trial_history(2);
    auto& qm = trial_history(3);
    auto& z = trial_history(4);

    const vec current_hfc(&current_history(5), bfc.size(), false, true);
    const vec current_hac(&current_history(5 + bfc.size()), bac.size(), false, true);
    const vec current_alpha(&current_history(5 + bfc.size() + bac.size()), ba.size(), false, true);
    const vec current_d(&current_history(5 + bfc.size() + bac.size() + ba.size()), bd.size(), false, true);
    vec hfc(&trial_history(5), bfc.size(), false, true);
    vec hac(&trial_history(5 + bfc.size()), bac.size(), false, true);
    vec alpha(&trial_history(5 + bfc.size() + bac.size()), ba.size(), false, true);
    vec d(&trial_history(5 + bfc.size() + bac.size() + ba.size()), bd.size(), false, true);

    iteration = 0.;
    const auto start_z = initial_check(current_z);

    // elastic unloading
    if(last_loading < -.5) return SUANPAN_SUCCESS;

    const auto ref_zr = trial_zr.mean();

    auto gamma = 0., ref_error = 0.;

    vec2 residual, incre;
    mat22 jacobian;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        auto km = 0., dkm = 0.;
        if(z > ref_zr) {
            if(const auto k_exp = std::exp((z - ref_zr) / (z - 1.) / kr); std::isfinite(k_exp)) {
                dkm = k_exp / kr * (1. - ref_zr) * std::pow(z - 1., -2.);
                km = 1. - k_exp;
            }
            else km = 1.;
        }

        q = current_q + gamma;
        qm = current_qm + km * gamma;

        const auto [u, du] = bound_u(qm, true);
        const auto pupg = du * km, pupz = du * gamma * dkm;

        const auto [hf, phfpg, phfpz] = compute_isotropic_bound(gamma, km, dkm);
        const auto [ha, phapg, phapz] = compute_kinematic_bound(gamma, km, dkm);

        vec top_alpha(ba.size(), fill::none), bot_alpha(ba.size(), fill::none);
        for(auto I = 0llu; I < ba.size(); ++I) {
            top_alpha(I) = ba[I].a() * gamma;
            bot_alpha(I) = 1. + ba[I].b() * gamma;
        }
        vec top_d(bd.size(), fill::none), bot_d(bd.size(), fill::none);
        for(auto I = 0llu; I < bd.size(); ++I) {
            top_d(I) = bd[I].a() * gamma;
            bot_d(I) = 1. + bd[I].b() * gamma;
        }

        const auto n = trial_stress(0) > ha * accu(current_alpha / bot_alpha) + (1. - z) * hf * accu(current_d / bot_d) ? 1. : -1.;

        const auto sum_alpha = accu(alpha = (top_alpha * n + current_alpha) / bot_alpha);
        const auto sum_d = accu(d = (top_d * n + current_d) / bot_d);

        auto palphapg = 0., pdpg = 0.;
        for(auto I = 0llu; I < ba.size(); ++I) palphapg += (ba[I].a() * n - ba[I].b() * alpha[I]) / bot_alpha[I];
        for(auto I = 0llu; I < bd.size(); ++I) pdpg += (bd[I].a() * n - bd[I].b() * d[I]) / bot_d[I];

        const auto trial_ratio = yield_ratio(z);
        const auto diff_z = z - start_z;

        residual(0) = std::fabs(trial_stress(0) - elastic * gamma * n - ha * sum_alpha + (z - 1.) * hf * sum_d) - z * hf;
        residual(1) = hf * diff_z - gamma * u * trial_ratio[0];

        jacobian(0, 0) = n * ((z - 1.) * (hf * pdpg + sum_d * phfpg) - ha * palphapg - sum_alpha * phapg) - elastic - z * phfpg;
        jacobian(0, 1) = n * ((hf + (z - 1.) * phfpz) * sum_d - phapz * sum_alpha) - hf - z * phfpz;

        jacobian(1, 0) = phfpg * diff_z - (u + gamma * pupg) * trial_ratio[0];
        jacobian(1, 1) = hf + phfpz * diff_z - gamma * (pupz * trial_ratio[0] + u * trial_ratio[1]);

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

        gamma = suanpan::clamp(gamma - incre(0), 0., 1.1 * abs_incre_strain);
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
