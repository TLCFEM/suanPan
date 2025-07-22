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

#include <Domain/DomainBase.h>
#include <Toolbox/tensor.h>

const double YLD0418P::root_two_third = std::sqrt(2. / 3.);
const span YLD0418P::sb{1, 6};
const mat YLD0418P::unit_dev_tensor = tensor::unit_deviatoric_tensor4v2();

YLD0418P::yield_t YLD0418P::compute_yield_surface(const vec3& psa, const mat33& pva, const vec3& psb, const mat33& pvb) const {
    auto f{0.};

    vec3 pfpa(fill::zeros), pfpb(fill::zeros), pfpaa(fill::zeros), pfpbb(fill::zeros);
    mat33 pfpab(fill::none);

    for(auto i = 0u; i < 3u; ++i)
        for(auto j = 0u; j < 3u; ++j) {
            const auto diff = (psa(i) - psb(j)) / ref_stress;
            const auto pow_diff = std::pow(std::max(datum::eps, std::fabs(diff)), exponent - 2.);
            f += pow_diff * diff * diff;
            pfpa(i) += pow_diff * diff;
            pfpb(j) -= pow_diff * diff;
            pfpaa(i) += pow_diff;
            pfpbb(j) += pow_diff;
            pfpab(i, j) = -pow_diff;
        }

    auto factor = exponent / ref_stress;

    pfpa *= factor;
    pfpb *= factor;

    factor *= (exponent - 1.) / ref_stress;

    pfpaa *= factor;
    pfpbb *= factor;
    pfpab *= factor;

    const mat66 proj_a = transform::eigen_to_tensor_base(pva).t() * C1;
    const mat66 proj_b = transform::eigen_to_tensor_base(pvb).t() * C2;

    const auto trans_a = proj_a.head_rows(3); // 3x6
    const auto trans_b = proj_b.head_rows(3); // 3x6

    const mat66 pfpmix = trans_a.t() * pfpab * trans_b;

    const auto kernel = [](const vec3& ps, const vec3& pfp) {
        const auto item = [&](const unsigned i, const unsigned j) { return 2. * (pfp(i) - pfp(j)) / (ps(i) - ps(j)); };

        return vec3{item(0, 1), item(1, 2), item(2, 0)};
    };

    return {f, trans_a.t() * pfpa + trans_b.t() * pfpb, proj_a.t() * diagmat(join_cols(pfpaa, kernel(psa, pfpa))) * proj_a + proj_b.t() * diagmat(join_cols(pfpbb, kernel(psb, pfpb))) * proj_b + pfpmix + pfpmix.t()};
}

YLD0418P::YLD0418P(const unsigned T, vec&& EE, vec&& VV, vec&& PP, const double M, const double RS, const unsigned HT, const double R)
    : DataYLD0418P{std::move(EE), std::move(VV), std::move(PP), M, std::fabs(RS)}
    , Material3D(T, R)
    , hardening_tag(HT) {
    C1.zeros();
    C2.zeros();

    C1(0, 1) = -parameter(0);
    C1(0, 2) = -parameter(1);
    C1(1, 0) = -parameter(2);
    C1(1, 2) = -parameter(3);
    C1(2, 0) = -parameter(4);
    C1(2, 1) = -parameter(5);

    C1(3, 3) = 2. * parameter(6);
    C1(4, 4) = 2. * parameter(7);
    C1(5, 5) = 2. * parameter(8);

    C2(0, 1) = -parameter(9);
    C2(0, 2) = -parameter(10);
    C2(1, 0) = -parameter(11);
    C2(1, 2) = -parameter(12);
    C2(2, 0) = -parameter(13);
    C2(2, 1) = -parameter(14);

    C2(3, 3) = 2. * parameter(15);
    C2(4, 4) = 2. * parameter(16);
    C2(5, 5) = 2. * parameter(17);
}

int YLD0418P::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(hardening_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", hardening_tag);
        return SUANPAN_FAIL;
    }

    hardening_expression = D->get_expression(hardening_tag);

    if(hardening_expression->input_size() != 1 || hardening_expression->output_size() != 1) {
        suanpan_error("The assigned expression {} is not a univariate function.\n", hardening_tag);
        return SUANPAN_FAIL;
    }

    if(!suanpan::approx_equal(1., hardening_expression->evaluate(0.).at(0))) {
        suanpan_error("The assigned expression {} does not evaluate to unity for trivial plastic strain.\n", hardening_tag);
        return SUANPAN_FAIL;
    }

    initial_stiffness = 1llu == modulus.n_elem && 1llu == ratio.n_elem ? tensor::isotropic_stiffness(modulus(0), ratio(0)) : tensor::orthotropic_stiffness(modulus, ratio);

    dev_ini_stiffness = unit_dev_tensor * (trial_stiffness = current_stiffness = initial_stiffness);

    initialize_history(7);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> YLD0418P::get_copy() { return std::make_unique<YLD0418P>(*this); }

int YLD0418P::update_trial_status(const vec& t_strain) {
    if(tensor::strain::norm(incre_strain = (trial_strain = t_strain) - current_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    const auto& current_ep = current_history(0);
    auto& ep = trial_history(0);
    vec plastic_strain(&trial_history(1), 6, false, true);

    const vec6 trial_dev_s = tensor::dev(trial_stress = (trial_stiffness = initial_stiffness) * (trial_strain - plastic_strain));

    auto gamma{0.};
    auto dev_s = trial_dev_s;

    vec7 incre(fill::none), residual(fill::none);
    mat77 jacobian(fill::none);

    auto counter{0u};
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        vec3 psa, psb;  // 3 eigenvalues
        mat33 pva, pvb; // 3x3 eigenvectors
        // !!! using strain version due to additional shear factors combined into transformation matrices
        if(!eig_sym(psa, pva, tensor::strain::to_tensor(C1 * dev_s), "std")) return SUANPAN_FAIL;
        // !!! using strain version due to additional shear factors combined into transformation matrices
        if(!eig_sym(psb, pvb, tensor::strain::to_tensor(C2 * dev_s), "std")) return SUANPAN_FAIL;

        const auto [f, pfps, pfpss] = compute_yield_surface(psa, pva, psb, pvb);
        const vec6 n = tensor::dev(pfps); // associated plastic flow direction
        const auto norm_n = root_two_third * tensor::strain::norm(n);
        const auto [k, dk] = compute_hardening(ep = current_ep + gamma * norm_n);

        residual(sa) = f - 4. * std::pow(k, exponent);

        if(1u == counter && residual(sa) < 0.) return SUANPAN_SUCCESS;

        const auto pk = -4. * exponent * std::pow(k, exponent - 1.) * dk;
        const vec6 en = initial_stiffness * n;

        residual(sb) = dev_s + en * gamma - trial_dev_s;

        jacobian(sa, sa) = pk * norm_n;
        jacobian(sa, sb) = pfps.t() + pk * gamma * two_third / norm_n * (n % tensor::strain::norm_weight).t() * pfpss;
        jacobian(sb, sa) = en;
        jacobian(sb, sb) = eye(6, 6) + gamma * dev_ini_stiffness * pfpss;

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        suanpan_debug("Local plasticity iteration error: {:.5E}.\n", error);
        if(error < tolerance * (1. + ref_stress) || (inf_norm(residual) < tolerance && counter > 5u)) {
            plastic_strain += n * gamma;

            trial_stress -= en * gamma;

            mat::fixed<7, 6> left(fill::none), right(fill::zeros);
            right.rows(sb) = dev_ini_stiffness;

            if(!solve(left, jacobian, right, solve_opts::equilibrate)) return SUANPAN_FAIL;

            trial_stiffness -= dev_ini_stiffness * join_rows(pfps, pfpss * gamma) * left;

            return SUANPAN_SUCCESS;
        }

        gamma = std::max(0., gamma - incre(sa));
        dev_s -= incre(sb);
    }
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

void YLD0418P::print() {
    if(1llu == modulus.n_elem && 1llu == ratio.n_elem) suanpan_info("YLD2004-18P model with E={:.5E} and nu={:.5E}.\n", modulus(0), ratio(0));
    else suanpan_info("YLD2004-18P model with E_1={:.5E}, E_2={:.5E}, E_3={:.5E}, G_{{12}}={:.5E}, G_{{23}}={:.5E}, G_{{13}}={:.5E}, and nu_{{12}}={:.5E}, nu_{{23}}={:.5E}, nu_{{13}}={:.5E}.\n", modulus(0), modulus(1), modulus(2), modulus(3), modulus(4), modulus(5), ratio(0), ratio(1), ratio(2));
}
