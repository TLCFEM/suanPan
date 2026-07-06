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

#include "Subloading1D.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

Subloading1D::Subloading1D(const unsigned T, DataSubloading1D&& D, const double R)
    : DataSubloading1D{std::move(D)}
    , Material1D(T, R) {}

int Subloading1D::initialize(const shared_ptr<DomainBase>& D) {
    if(nullptr != D) incre_time = &D->get_factory()->modify_incre_time();

    trial_stiffness = current_stiffness = initial_stiffness = elastic;

    initialize_history(4u + static_cast<unsigned>(b.size() + c.size()));

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Subloading1D::unique_copy() { return std::make_unique<Subloading1D>(*this); }

int Subloading1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    const auto abs_incre_strain = std::fabs(incre_strain(0));
    const auto zero_increment = abs_incre_strain <= datum::eps;

    if(!is_viscous && zero_increment) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    const auto& current_q = current_history(1);
    const auto& current_z = current_history(2);
    auto& iteration = trial_history(0);
    auto& q = trial_history(1);
    auto& z = trial_history(2);
    auto& zv = trial_history(3);

    const auto current_na = b.empty() ? vec{} : vec(&current_history(4), static_cast<uword>(b.size()), false, true);
    const auto current_nd = c.empty() ? vec{} : vec(&current_history(4 + static_cast<uword>(b.size())), static_cast<uword>(c.size()), false, true);
    auto na = b.empty() ? vec{} : vec(&trial_history(4), static_cast<uword>(b.size()), false, true);
    auto nd = c.empty() ? vec{} : vec(&trial_history(4 + static_cast<uword>(b.size())), static_cast<uword>(c.size()), false, true);

    const auto norm_mu = mu / (incre_time && *incre_time > 0. ? *incre_time : 1.);

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

        const auto [y, dy] = iso_bound(q, true);
        const auto [a, da] = kin_bound(q, true);

        vec bot_na(static_cast<uword>(b.size()), fill::none), bot_nd(static_cast<uword>(c.size()), fill::none);
        for(uword I{0}; I < b.size(); ++I) bot_na(I) = 1. + b[I].r() * gamma;
        for(uword I{0}; I < c.size(); ++I) bot_nd(I) = 1. + c[I].r() * gamma;

        const auto n = trial_stress(0) - a * accu(current_na / bot_na) + (z - 1.) * y * accu(current_nd / bot_nd) > 0. ? 1. : -1.;

        for(uword I{0}; I < b.size(); ++I) na(I) = (b[I].rb() * gamma * n + current_na(I)) / bot_na(I);
        for(uword I{0}; I < c.size(); ++I) nd(I) = (c[I].rb() * gamma * n + current_nd(I)) / bot_nd(I);

        const auto sum_na = accu(na), sum_nd = accu(nd);

        if(1u == counter && !zero_increment) {
            const auto s = (y * sum_nd + a * sum_na - current_stress(0)) / (trial_stress(0) - current_stress(0));
            if(s >= 1.) {
                zv = ((trial_stress(0) - a * sum_na) / y - sum_nd) / (n - sum_nd);
                if(zv < z) {
                    // elastic unloading
                    z = zv;
                    return SUANPAN_SUCCESS;
                }
                // return SUANPAN_SUCCESS;
            }
            if(s > 0.) start_z = 0.;
        }

        auto dna = 0., dnd = 0.;
        for(uword I{0}; I < b.size(); ++I) dna += b[I].r() * (b[I].b() * n - na[I]) / bot_na[I];
        for(uword I{0}; I < c.size(); ++I) dnd += c[I].r() * (c[I].b() * n - nd[I]) / bot_nd[I];

        const auto trial_ratio = yield_ratio(z);
        const auto avg_rate = u * trial_ratio[0];
        const auto fraction_term = (cv * z - zv) * norm_mu * gamma + 1.;
        const auto power_term = std::pow(fraction_term, nv - 1.);

        residual(0) = std::fabs(trial_stress(0) - elastic * gamma * n - a * sum_na + (z - 1.) * y * sum_nd) - zv * y;
        residual(1) = z - start_z - gamma * avg_rate;
        residual(2) = zv - fraction_term * power_term * z;

        jacobian(0, 0) = n * ((z - 1.) * (y * dnd + sum_nd * dy) - (a * dna + sum_na * da)) - elastic - zv * dy;
        jacobian(0, 1) = -y;
        jacobian(0, 2) = n * y * sum_nd;

        jacobian(1, 0) = -avg_rate;
        jacobian(1, 1) = 0.;
        jacobian(1, 2) = 1. - u * gamma * trial_ratio[1];

        jacobian(2, 0) = -z * nv * power_term * (cv * z - zv) * norm_mu;
        jacobian(2, 1) = 1. + z * nv * power_term * norm_mu * gamma;
        jacobian(2, 2) = -power_term * (fraction_term + z * nv * cv * norm_mu * gamma);

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

        const auto error = suanpan::inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || suanpan::inf_norm(residual) < tolerance) && counter > 5u)) {
            iteration = counter;
            trial_stress -= elastic * gamma * n;
            trial_stiffness += elastic / det(jacobian) * elastic * det(jacobian.submat(1, 1, 2, 2));
            return SUANPAN_SUCCESS;
        }

        gamma -= incre(0);
        if(gamma < 0.) gamma = 0.;

        z = suanpan::clamp_unit(z - incre(2));
        zv = is_viscous ? suanpan::clamp(zv - incre(1), z, cv * z) : z;
    }
}

int Subloading1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int Subloading1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int Subloading1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void Subloading1D::print() {
    suanpan_info("A uniaxial combined hardening material using subloading surface model with optional viscosity. doi:10.1007/s00707-025-04339-0\n");
    Material1D::print();
}
