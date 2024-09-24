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

#include "Subloading1D.h"

const double Subloading1D::rate_bound = log(r_bound);

Subloading1D::Subloading1D(const unsigned T, DataSubloading1D&& D, const double R)
    : DataSubloading1D{std::move(D)}
    , Material1D(T, R) {}

int Subloading1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic;

    initialize_history(3);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Subloading1D::get_copy() { return make_unique<Subloading1D>(*this); }

int Subloading1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    const auto& current_r = current_history(2);
    auto& back_stress = trial_history(0);
    auto& plastic_strain = trial_history(1);
    auto& r = trial_history(2);

    const auto eta = trial_stress(0) - back_stress;
    const auto norm_eta = fabs(eta);

    auto gamma = 0., ref_error = 0.;

    vec2 residual, incre;
    mat22 jacobian;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto iso_stress = yield + isotropic * (plastic_strain + gamma);
        const auto r_rate = -u * (r > r_bound ? log(r) : rate_bound);
        residual(0) = norm_eta - (elastic + kinematic) * gamma - r * iso_stress;
        residual(1) = r - current_r - r_rate * gamma;

        jacobian(0, 0) = -elastic - kinematic - r * isotropic;
        jacobian(0, 1) = -iso_stress;
        jacobian(1, 0) = -r_rate;
        jacobian(1, 1) = 1. + gamma * (r > r_bound ? u / r : 0.);

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) {
            if(gamma > 0.) {
                plastic_strain += gamma;
                if(eta < 0.) gamma = -gamma;
                back_stress += kinematic * gamma;
                trial_stress -= elastic * gamma;

                const vec2 right = solve(jacobian, vec2{elastic, 0.});
                trial_stiffness += elastic * right(0);
            }
            else r = norm_eta / iso_stress;
            break;
        }

        gamma -= incre(0);
        r -= incre(1);
    }

    return SUANPAN_SUCCESS;
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
    suanpan_info("A uniaxial combined hardening material using subloading surface model.\n");
    Material1D::print();
}
