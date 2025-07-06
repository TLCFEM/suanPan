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

#include "YLD0418P.h"

#include <Toolbox/ridders.hpp>
#include <Toolbox/tensor.h>

const double YLD0418P::root_two_third = std::sqrt(2. / 3.);
const span YLD0418P::sb{1, 6};
const mat YLD0418P::unit_dev_tensor = tensor::unit_deviatoric_tensor4v2();

YLD0418P::YLD0418P(const unsigned T, vec&& EE, vec&& VV, vec&& PP, const double M, const double R)
    : DataYLD0418P{std::move(EE), std::move(VV), std::move(PP), M}
    , Material3D(T, R) {
    C1.zeros();
    C2.zeros();

    C1(0, 1) = -parameter(0);
    C1(0, 2) = -parameter(1);
    C1(1, 0) = -parameter(2);
    C1(1, 2) = -parameter(3);
    C1(2, 0) = -parameter(4);
    C1(2, 1) = -parameter(5);

    C1(3, 3) = parameter(6);
    C1(4, 4) = parameter(7);
    C1(5, 5) = parameter(8);

    C2(0, 1) = -parameter(9);
    C2(0, 2) = -parameter(10);
    C2(1, 0) = -parameter(11);
    C2(1, 2) = -parameter(12);
    C2(2, 0) = -parameter(13);
    C2(2, 1) = -parameter(14);

    C2(3, 3) = parameter(15);
    C2(4, 4) = parameter(16);
    C2(5, 5) = parameter(17);

    C1 *= unit_dev_tensor;
    C2 *= unit_dev_tensor;
}

int YLD0418P::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::orthotropic_stiffness(modulus, ratio);

    initialize_history(7);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> YLD0418P::get_copy() { return std::make_unique<YLD0418P>(*this); }

int YLD0418P::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    const auto norm_incre_strain = tensor::strain::norm(incre_strain);

    if(norm_incre_strain <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    const auto& current_ep = current_history(0);
    auto& ep = trial_history(0);
    vec plastic_strain(&trial_history(1), 6, false, true);

    const vec6 predictor_s = trial_stress = (trial_stiffness = initial_stiffness) * (trial_strain - plastic_strain);

    auto gamma{0.};

    vec7 incre(fill::none), residual(fill::none);
    mat77 jacobian(fill::none);

    auto counter{0u};
    auto try_bisection{false};
    while(true) {
        if(max_iteration == ++counter) {
            try_bisection = true;
            return SUANPAN_FAIL;
        }

        vec3 principal_a, principal_b;    // 3
        mat33 principal_da, principal_db; // 3x3
        if(!eig_sym(principal_a, principal_da, tensor::stress::to_tensor(C1 * trial_stress), "std")) return SUANPAN_FAIL;
        if(!eig_sym(principal_b, principal_db, tensor::stress::to_tensor(C2 * trial_stress), "std")) return SUANPAN_FAIL;

        const mat trans_ra = transform::compute_jacobian_nominal_to_principal(principal_da) * C1; // 3x6
        const mat trans_rb = transform::compute_jacobian_nominal_to_principal(principal_db) * C2; // 3x6

        auto f{0.};

        rowvec3 dfdpa(fill::zeros), dfdpb(fill::zeros);
        vec3 dfdaa(fill::zeros), dfdbb(fill::zeros);
        mat33 dfdab(fill::zeros);
        for(auto i = 0u; i < 3u; ++i) {
            for(auto j = 0u; j < 3u; ++j) {
                const auto diff = (principal_a(i) - principal_b(j)) / ref_stress;
                const auto abs_diff = std::fabs(diff);
                const auto pow_diff = std::pow(std::max(datum::eps, abs_diff), exponent - 2.);
                f += pow_diff * abs_diff * abs_diff;
                dfdpa(i) += pow_diff * diff;
                dfdpb(j) -= pow_diff * diff;
                dfdaa(i) += pow_diff;
                dfdbb(j) += pow_diff;
                dfdab(i, j) -= pow_diff;
            }
        }
        dfdpa *= exponent / ref_stress;
        dfdpb *= exponent / ref_stress;
        dfdaa *= exponent * (exponent - 1.) / ref_stress / ref_stress;
        dfdbb *= exponent * (exponent - 1.) / ref_stress / ref_stress;
        dfdab *= exponent * (exponent - 1.) / ref_stress / ref_stress;

        const vec6 dfds = (dfdpa * trans_ra + dfdpb * trans_rb).t();
        mat66 dft = trans_ra.t() * dfdab * trans_rb;
        mat66 dfdss = trans_ra.t() * diagmat(dfdaa) * trans_ra + trans_rb.t() * diagmat(dfdbb) * trans_rb + dft + dft.t();

        const auto norm_n = root_two_third * tensor::strain::norm(dfds);

        const auto k = compute_k(ep = current_ep + gamma * norm_n);
        const auto pk = -4. * exponent * std::pow(k, exponent - 1.) * compute_dk(ep);

        residual(sa) = f - 4. * std::pow(k, exponent);
        residual(sb) = trial_stress + gamma * initial_stiffness * dfds - predictor_s;

        jacobian(sa, sa) = pk * norm_n;
        jacobian(sa, sb) = dfds.t() + pk * gamma * two_third / norm_n * (dfds % tensor::strain::norm_weight).t() * dfdss;
        jacobian(sb, sa) = initial_stiffness * dfds;
        jacobian(sb, sb) = eye(6, 6) + gamma * initial_stiffness * dfdss;

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        suanpan_debug("Local plasticity iteration error: {:.5E}.\n", error);
        if(error < tolerance * (1. + ref_stress) || (inf_norm(residual) < tolerance && counter > 5u) || try_bisection) {
            plastic_strain += gamma * dfds;

            mat::fixed<7, 6> left(fill::none), right(fill::zeros);
            right.rows(sb) = initial_stiffness;

            if(!solve(left, jacobian, right, solve_opts::equilibrate)) return SUANPAN_FAIL;

            trial_stiffness = left.rows(sb);

            return SUANPAN_SUCCESS;
        }

        gamma = std::max(0., gamma - incre(sa));
        trial_stress -= incre(sb);
    }

    return SUANPAN_SUCCESS;
}

int YLD0418P::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int YLD0418P::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int YLD0418P::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void YLD0418P::print() { suanpan_info("YLD2004-18P model with E_1={:.5E}, E_2={:.5E}, E_3={:.5E}, G_{{12}}={:.5E}, G_{{23}}={:.5E}, G_{{13}}={:.5E}, and nu_{{12}}={:.5E}, nu_{{23}}={:.5E}, nu_{{13}}={:.5E}.\n", modulus(0), modulus(1), modulus(2), modulus(3), modulus(4), modulus(5), ratio(0), ratio(1), ratio(2)); }
