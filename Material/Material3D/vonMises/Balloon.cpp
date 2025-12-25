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

#include "Balloon.h"

#include <Toolbox/tensor.h>

const mat Balloon::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

auto Balloon ::compute_isotropic_bound(const double incre_q, const double km, const double dkm) {
    const auto& qm = trial_history(3);

    const auto current_hfc = bfc.empty() ? vec{} : vec(&current_history(5), bfc.size(), false, true);
    auto hfc = bfc.empty() ? vec{} : vec(&trial_history(5), bfc.size(), false, true);

    const auto kc = 1. - km, dkc = -dkm;

    const auto [fm, dfm] = bound_fm(qm, true);
    const auto [fc, dfc] = bound_fc(qm, false);

    auto phfpg = dfm * root_two_third * km, phfpz = dfm * incre_q * dkm;
    const auto pfcpg = dfc * root_two_third * km, pfcpz = dfc * incre_q * dkm;
    const auto incre_qc = kc * incre_q;
    for(auto I = 0llu; I < bfc.size(); ++I) {
        const auto bot_fc = 1. + bfc[I].b() * incre_qc;
        hfc(I) = (bfc[I].a() * incre_qc * fc + current_hfc(I)) / bot_fc;
        phfpg += (bfc[I].a() * (kc * root_two_third * fc + incre_qc * pfcpg) - hfc(I) * bfc[I].b() * kc * root_two_third) / bot_fc;
        phfpz += (bfc[I].a() * (dkc * incre_q * fc + incre_qc * pfcpz) - hfc(I) * bfc[I].b() * dkc * incre_q) / bot_fc;
    }

    if(const auto hf = fm + accu(hfc); hf > 0.) return std::make_tuple(hf, phfpg, phfpz);

    return std::make_tuple(0., 0., 0.);
}

auto Balloon ::compute_kinematic_bound(const double incre_q, const double km, const double dkm) {
    const auto& qm = trial_history(3);

    const auto current_hac = bac.empty() ? vec{} : vec(&current_history(5 + bfc.size()), bac.size(), false, true);
    auto hac = bac.empty() ? vec{} : vec(&trial_history(5 + bfc.size()), bac.size(), false, true);

    const auto kc = 1. - km, dkc = -dkm;

    const auto [am, dam] = bound_am(qm, true);
    const auto [ac, dac] = bound_ac(qm, false);

    auto phapg = dam * root_two_third * km, phapz = dam * incre_q * dkm;
    const auto pacpg = dac * root_two_third * km, pacpz = dac * incre_q * dkm;
    const auto incre_qc = kc * incre_q;
    for(auto I = 0llu; I < bac.size(); ++I) {
        const auto bot_ac = 1. + bac[I].b() * incre_qc;
        hac(I) = (bac[I].a() * incre_qc * ac + current_hac(I)) / bot_ac;
        phapg += (bac[I].a() * (kc * root_two_third * ac + incre_qc * pacpg) - hac(I) * bac[I].b() * kc * root_two_third) / bot_ac;
        phapz += (bac[I].a() * (dkc * incre_q * ac + incre_qc * pacpz) - hac(I) * bac[I].b() * dkc * incre_q) / bot_ac;
    }

    if(const auto ha = am + accu(hac); ha > 0.) return std::make_tuple(ha, phapg, phapz);

    return std::make_tuple(0., 0., 0.);
}

Balloon::Balloon(const unsigned T, DataBalloon&& D, const double R)
    : DataBalloon{std::move(D)}
    , Material3D(T, R) { access::rw(tolerance) = 1E-13; }

int Balloon::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic, poisson);

    initialize_history(5u + static_cast<unsigned>(bfc.size() + bac.size() + 6u * (bna.size() + bnd.size())));

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Balloon::get_copy() { return std::make_unique<Balloon>(*this); }

int Balloon::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    const auto norm_incre_strain = tensor::strain::norm(incre_strain);
    if(norm_incre_strain <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    const auto& current_q = current_history(2);
    const auto& current_qm = current_history(3);
    const auto& current_z = current_history(4);
    auto& iteration = trial_history(0);
    auto& last_loading = trial_history(1);
    auto& q = trial_history(2);
    auto& qm = trial_history(3);
    auto& z = trial_history(4);

    const auto offset_na = 5u + bfc.size() + bac.size();
    const auto offset_nd = offset_na + 6u * bna.size();

    const vec trial_s = tensor::dev(trial_stress);

    auto ref_zr = trial_zr(zr_type);

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

        auto km = 0., dkm = 0.;
        if(z > ref_zr) {
            if(const auto k_exp = std::exp((z - ref_zr) / (z - 1.) / kr); std::isfinite(k_exp)) {
                dkm = k_exp / kr * (1. - ref_zr) * std::pow(z - 1., -2.);
                km = 1. - k_exp;
            }
            else km = 1.;
        }

        const auto incre_q = root_two_third * gamma;
        q = current_q + incre_q;
        qm = current_qm + km * incre_q;

        const auto [u, du] = bound_u(qm, true);
        const auto pupg = du * root_two_third * km, pupz = du * incre_q * dkm;

        const auto [hf, phfpg, phfpz] = compute_isotropic_bound(incre_q, km, dkm);
        const auto [ha, phapg, phapz] = compute_kinematic_bound(incre_q, km, dkm);

        vec6 sum_na(fill::zeros);
        vec top_na(bna.size(), fill::none), bot_na(bna.size(), fill::none);
        auto dna{0.};
        for(auto I = 0llu; I < bna.size(); ++I) {
            top_na(I) = bna[I].a() * gamma;
            sum_na += vec(&current_history(offset_na + 6u * I), 6, false, true) / (bot_na(I) = 1. + bna[I].b() * gamma);
            dna += (bna[I].a() - top_na(I) / bot_na(I) * bna[I].b()) / bot_na(I);
        }
        vec6 sum_nd(fill::zeros);
        vec top_nd(bnd.size(), fill::none), bot_nd(bnd.size(), fill::none);
        auto dnd{0.};
        for(auto I = 0llu; I < bnd.size(); ++I) {
            top_nd(I) = bnd[I].a() * gamma;
            sum_nd += vec(&current_history(offset_nd + 6u * I), 6, false, true) / (bot_nd(I) = 1. + bnd[I].b() * gamma);
            dnd += (bnd[I].a() - top_nd(I) / bot_nd(I) * bnd[I].b()) / bot_nd(I);
        }

        if(1u == counter) {
            const vec incre_s = trial_s - tensor::dev(current_stress);
            const vec ref = trial_s - ha * sum_na - hf * sum_nd;
            const vec base = ref - incre_s;

            const auto aa = two_third - tensor::stress::double_contraction(sum_nd);
            const auto bb = tensor::stress::double_contraction(sum_nd, ref);
            const auto cc = tensor::stress::double_contraction(ref);

            const auto incre_incre = tensor::stress::double_contraction(incre_s);
            const auto incre_d = tensor::stress::double_contraction(incre_s, sum_nd);

            auto x = .5;
            auto inner_counter = 0u;
            while(true) {
                if(max_iteration == ++inner_counter) {
                    suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
                    return SUANPAN_FAIL;
                }

                const vec middle = base + x * incre_s;
                const auto middle_d = tensor::stress::double_contraction(sum_nd, middle);
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
                    if(std::signbit(last_loading) == std::signbit(x)) ref_zr = trial_zr.enqueue(z)(zr_type);
                    if(x >= 1.) {
                        // elastic unloading
                        last_loading = -1.;
                        z = (bb + std::sqrt(bb * bb + aa * cc)) / aa / hf;
                        return SUANPAN_SUCCESS;
                    }
                    last_loading = 1.;
                    if(x > 0.) {
                        start_z = (middle_d + tmp_sqrt) / aa / hf;
                        suanpan_debug("Initial yield ratio: {:.5E}, corrected yield ratio: {:.5E}.\n", current_z, start_z);
                    }
                    break;
                }

                x -= incre_x;
            }
        }

        const vec zeta = trial_s - ha * sum_na + (z - 1.) * hf * sum_nd;
        const auto norm_zeta = tensor::stress::norm(zeta);
        const vec n = zeta / norm_zeta;

        // update history variables
        // they are not used in the state determination algorithm for the current step (but next one)
        for(auto I = 0llu; I < bna.size(); ++I) vec(&trial_history(offset_na + 6u * I), 6, false, true) = (top_na(I) * n + vec(&current_history(offset_na + 6u * I), 6, false, true)) / bot_na(I);
        for(auto I = 0llu; I < bnd.size(); ++I) vec(&trial_history(offset_nd + 6u * I), 6, false, true) = (top_nd(I) * n + vec(&current_history(offset_nd + 6u * I), 6, false, true)) / bot_nd(I);

        const vec pzetapz = hf * sum_nd + (z - 1.) * phfpz * sum_nd - phapz * sum_na;
        const vec pzetapg = (z - 1.) * phfpg * sum_nd - phapg * sum_na; // just a part of it

        sum_na.zeros();
        sum_nd.zeros();
        for(auto I = 0llu; I < bna.size(); ++I) sum_na -= bna[I].b() * std::pow(bot_na(I), -2.) * vec(&current_history(offset_na + 6u * I), 6, false, true);
        for(auto I = 0llu; I < bnd.size(); ++I) sum_nd -= bnd[I].b() * std::pow(bot_nd(I), -2.) * vec(&current_history(offset_nd + 6u * I), 6, false, true);

        const auto trial_ratio = yield_ratio(z);
        const auto diff_z = z - start_z;

        const auto factor_na = accu(top_na / bot_na), factor_nd = accu(top_nd / bot_nd);
        residual(0) = norm_zeta - double_shear * gamma - ha * factor_na + (z - 1.) * hf * factor_nd - root_two_third * hf * z;
        residual(1) = hf * diff_z - gamma * u * trial_ratio[0];

        jacobian(0, 0) = tensor::stress::double_contraction(n, pzetapg + (z - 1.) * hf * sum_nd - ha * sum_na) - double_shear - phapg * factor_na - ha * dna + (z - 1.) * (phfpg * factor_nd + hf * dnd) - root_two_third * phfpg * z;
        jacobian(0, 1) = tensor::stress::double_contraction(n, pzetapz) - phapz * factor_na + (hf + (z - 1.) * phfpz) * factor_nd - root_two_third * (hf + phfpz * z);
        jacobian(1, 0) = phfpg * diff_z - (u + gamma * pupg) * trial_ratio[0];
        jacobian(1, 1) = hf + phfpz * diff_z - gamma * (pupz * trial_ratio[0] + u * trial_ratio[1]);

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

int Balloon::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int Balloon::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int Balloon::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void Balloon::print() {
    suanpan_info("The 3D version of the Balloon-v1 model.\n");
}
