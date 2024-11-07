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

#include "VAFCRP1D.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

std::tuple<double, double> VAFCRP1D::isotropic_bound(const double p) const {
    const auto exp_term = saturated * exp(-m * p);

    auto k = yield + saturated + hardening * p - exp_term;
    auto dk = hardening + m * exp_term;
    if(k < 0.) k = dk = 0.;
    return {k, dk};
}

bool VAFCRP1D::is_elastic(const double stress, const vec& history) const {
    const auto [k, dk] = isotropic_bound(history(size));
    return fabs(stress - accu(history.head(size))) < k;
}

int VAFCRP1D::partial_loading(double& fragment, const vec& start_history, const double start_stress, const double diff_stress) {
    const auto& start_p = start_history(size);
    const vec start_beta(&start_history(0), size);

    auto& p = trial_history(size);
    vec beta(&trial_history(0), size, false, true);

    auto inter_history = start_history;
    auto& inter_p = inter_history(size);
    vec inter_beta(&inter_history(0), size, false, true);

    const auto norm_mu = mu / (incre_time && *incre_time > 0. ? *incre_time : 1.);

    auto gamma = 0., ref_error = 1.;

    vec2 residual, incre;
    mat22 jacobian;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        residual.zeros();
        jacobian.zeros();

        if(fragment > 1.) fragment = 1.;

        const auto inter_gamma = middle_point * gamma;

        p = start_p + gamma;
        inter_p = start_p + inter_gamma;

        const auto [inter_k, inter_dk] = isotropic_bound(inter_p);

        const vec inter_bottom = 1. + b * inter_gamma;

        const auto inter_trial_stress = start_stress + middle_point * fragment * diff_stress;
        const auto inter_n = inter_trial_stress > accu(start_beta / inter_bottom) ? 1. : -1.;

        inter_beta = (start_beta + inter_n * inter_gamma * a) / inter_bottom;
        const auto inter_dbeta = accu((a * inter_n - b % start_beta) / square(inter_bottom));

        beta = (inter_beta - (1. - middle_point) * start_beta) / middle_point;

        const auto inter_stress = inter_trial_stress - inter_n * elastic_modulus * inter_gamma;
        const auto inter_eta = inter_stress - accu(inter_beta);

        const auto fraction_term = norm_mu * gamma + 1.;
        const auto power_term = pow(fraction_term, epsilon - 1.);

        residual(0) = fabs(inter_eta) - fraction_term * power_term * inter_k;

        jacobian(0, 0) = -middle_point * (elastic_modulus + inter_n * inter_dbeta) - power_term * (norm_mu * epsilon * inter_k + fraction_term * inter_dk * middle_point);

        if(fragment >= 1.) jacobian(1, 1) = 1.;
        else {
            const auto eta = start_stress + fragment * diff_stress - inter_n * elastic_modulus * gamma - accu(beta);
            const auto n = eta > 0. ? 1. : -1.;

            const auto [k, dk] = isotropic_bound(p);

            residual(1) = fabs(eta) - k;

            jacobian(0, 1) = inter_n * middle_point * diff_stress;

            jacobian(1, 0) = -n * (inter_n * elastic_modulus + inter_dbeta) - dk;
            jacobian(1, 1) = n * diff_stress;
        }

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance * yield) && counter > 5u)) {
            trial_stress = start_stress + fragment * diff_stress - inter_n * elastic_modulus * gamma;

            if(fragment >= 1.) trial_stiffness = elastic_modulus + elastic_modulus / jacobian(0, 0) * elastic_modulus;

            return SUANPAN_SUCCESS;
        }

        gamma -= incre(0);
        fragment -= incre(1);
    }
}

int VAFCRP1D::partial_loading(const vec& start_history, const double start_stress) {
    if(is_elastic(start_stress, start_history)) return SUANPAN_SUCCESS;

    trial_history = start_history;
    auto& p = trial_history(size);
    vec beta(&trial_history(0), size, false, true);

    const auto norm_mu = mu / (incre_time && *incre_time > 0. ? *incre_time : 1.);

    auto gamma = 0., ref_error = 1.;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto [k, dk] = isotropic_bound(p);

        const vec bottom = 1. + b * gamma;

        const auto xi = start_stress - accu(beta / bottom);
        const auto q = fabs(xi) - (elastic_modulus + accu(a / bottom)) * gamma;

        const auto fraction_term = norm_mu * gamma + 1.;
        const auto power_term = pow(fraction_term, epsilon - 1.);

        auto jacobian = -elastic_modulus - power_term * (norm_mu * epsilon * k + fraction_term * dk);
        if(xi > 0.) jacobian += accu((b % beta - a) / square(bottom));
        else jacobian -= accu((b % beta + a) / square(bottom));

        const auto residual = q - fraction_term * power_term * k;

        const auto incre = residual / jacobian;
        const auto error = fabs(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || fabs(residual) < tolerance * yield) && counter > 5u)) {
            if(xi > 0.) {
                beta += a * gamma;

                trial_stress = start_stress - elastic_modulus * gamma;
            }
            else {
                beta -= a * gamma;

                trial_stress = start_stress + elastic_modulus * gamma;
            }

            beta /= bottom;

            trial_stiffness = elastic_modulus + elastic_modulus / jacobian * elastic_modulus;

            return SUANPAN_SUCCESS;
        }

        gamma -= incre;
        p -= incre;
    }
}

VAFCRP1D::VAFCRP1D(const unsigned T, DataVAFCRP1D&& D, const double R)
    : DataVAFCRP1D(std::move(D))
    , Material1D(T, R) { access::rw(tolerance) = 1E-15; }

int VAFCRP1D::initialize(const shared_ptr<DomainBase>& D) {
    if(nullptr != D) incre_time = &D->get_factory()->modify_incre_time();

    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(1 + size);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> VAFCRP1D::get_copy() { return make_unique<VAFCRP1D>(*this); }

double VAFCRP1D::get_parameter(const ParameterType P) const {
    if(ParameterType::ELASTICMODULUS == P) return elastic_modulus;
    return 0.;
}

int VAFCRP1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;
    incre_stress = (trial_stiffness = initial_stiffness) * incre_strain;
    trial_stress = current_stress + incre_stress;

    // pure loading
    if(const auto current_eta = current_stress(0) - accu(current_history.head(size)); current_eta * incre_stress(0) >= 0. || is_elastic(current_stress(0), current_history)) return partial_loading(current_history, trial_stress(0));

    auto fragment = 0.;
    if(SUANPAN_SUCCESS != partial_loading(fragment, current_history, current_stress(0), incre_stress(0))) return SUANPAN_FAIL;
    if(fragment >= 1.) return SUANPAN_SUCCESS;

    trial_stress += (1. - fragment) * incre_stress;

    return partial_loading(trial_history, trial_stress(0));
}

int VAFCRP1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int VAFCRP1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int VAFCRP1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void VAFCRP1D::print() {
    suanpan_info("A uniaxial VAFCRP material model.\n");
    Material1D::print();
}
