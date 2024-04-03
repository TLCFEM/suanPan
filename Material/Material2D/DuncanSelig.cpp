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

#include "DuncanSelig.h"
#include <Domain/DomainBase.h>

double DuncanSelig::dev(const vec& t_stress) { return std::sqrt(std::pow(t_stress(0) - t_stress(1), 2) + std::pow(2. * t_stress(2), 2)); }

rowvec3 DuncanSelig::der_dev(const vec& t_stress) { return {t_stress(0) - t_stress(1), t_stress(1) - t_stress(0), 4. * t_stress(2)}; }

int DuncanSelig::project_to_surface(double& elastic_portion) {
    const auto max_dev_stress = trial_history(0);

    const vec elastic_stress = initial_stiffness * incre_strain;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const vec t_stress = current_stress + elastic_stress * elastic_portion;
        const auto dev_stress = dev(t_stress);
        const auto residual = dev_stress - max_dev_stress;
        const auto incre = residual * dev_stress / dot(der_dev(t_stress), elastic_stress);

        if(std::abs(incre) < tolerance || (std::abs(residual) < tolerance && counter > 5u)) break;

        elastic_portion -= incre;
    }

    return SUANPAN_SUCCESS;
}

std::tuple<double, double, rowvec3, rowvec3> DuncanSelig::compute_moduli() {
    // principal stresses

    const auto center = -.5 * (trial_stress(0) + trial_stress(1));
    const rowvec3 dcds{-.5, -.5, 0.};

    auto radius = .5 * dev(trial_stress);
    rowvec3 drds(fill::zeros);
    if(radius > datum::eps) drds = der_dev(trial_stress) / radius * .25;

    const auto s1 = center + radius, s3 = center - radius;

    const rowvec3 ds1ds = dcds + drds, ds3ds = dcds - drds;

    // for elastic modulus
    double phi, dphids3;
    if(s3 < p_atm) {
        phi = ini_phi;
        dphids3 = 0.;
    }
    else {
        phi = ini_phi - ten_fold_phi_diff * log10(s3 / p_atm);
        dphids3 = -ten_fold_phi_diff / (s3 * log(10));
        if(phi < 0.) phi = dphids3 = 0.;
    }

    static constexpr auto min_ratio = 1.;

    const auto denom = 1. - std::sin(phi);
    auto max_dev_stress = 2. / r_f * (cohesion * std::cos(phi) + s3 * std::sin(phi)) / denom;
    auto dmdsds3 = 0.;
    if(max_dev_stress > min_ratio * p_atm) {
        const auto pmdspphi = 2. / r_f * (s3 * std::cos(phi) / denom / denom + cohesion / denom);
        const auto pmdsps3 = 2. / r_f * std::sin(phi) / denom;
        dmdsds3 = pmdspphi * dphids3 + pmdsps3;
    }
    else max_dev_stress = min_ratio * p_atm;

    const auto dev_stress = s1 - s3;
    const auto pdsps1 = 1.;
    const auto pdsps3 = -1.;

    double ini_elastic, deids3;
    if(s3 < -min_ratio * p_atm) {
        ini_elastic = ref_elastic * std::pow(.01, n);
        deids3 = 0.;
    }
    else if(s3 < min_ratio * p_atm) {
        ini_elastic = ref_elastic * std::pow(min_ratio, n);
        deids3 = 0.;
    }
    else {
        ini_elastic = ref_elastic * std::pow(s3 / p_atm, n);
        deids3 = n * ini_elastic / s3;
    }

    const auto pepei = std::pow(1. - dev_stress / max_dev_stress, 2.);
    const auto elastic = ini_elastic * pepei;
    const auto pepds = -2. * ini_elastic * (1. - dev_stress / max_dev_stress) / max_dev_stress;
    const auto pepmds = 2. * ini_elastic * (1. - dev_stress / max_dev_stress) * dev_stress / max_dev_stress / max_dev_stress;

    const auto peps1 = pepds * pdsps1;
    const auto peps3 = pepei * deids3 + pepds * pdsps3 + pepmds * dmdsds3;

    const rowvec3 deds = peps1 * ds1ds + peps3 * ds3ds;

    // for bulk modulus

    double bulk, pkps3;
    if(s3 < -min_ratio * p_atm) {
        bulk = ref_bulk * std::pow(.01, m);
        pkps3 = 0.;
    }
    else if(s3 < min_ratio * p_atm) {
        bulk = ref_bulk * std::pow(min_ratio, m);
        pkps3 = 0.;
    }
    else {
        bulk = ref_bulk * std::pow(s3 / p_atm, m);
        pkps3 = m * bulk / s3;
    }

    rowvec3 dkds = pkps3 * ds3ds;

    if(3. * bulk < elastic) {
        bulk = elastic / 3.;
        dkds = deds / 3.;
    }

    return {elastic, bulk, deds, dkds};
}

DuncanSelig::DuncanSelig(const unsigned T, const vec& P, const double R)
    : Material2D(T, PlaneType::E, R)
    , p_atm(P(0))
    , ref_elastic(P(1))
    , n(P(2))
    , ref_bulk(P(3))
    , m(P(4))
    , ini_phi(P(5))
    , ten_fold_phi_diff(P(6))
    , r_f(P(7))
    , cohesion(P(8)) { access::rw(tolerance) = 1E-13; }

int DuncanSelig::initialize(const shared_ptr<DomainBase>&) {
    const auto [elastic, bulk, deds, dkds] = compute_moduli();

    initial_stiffness.zeros(3, 3);
    initial_stiffness(0, 0) = initial_stiffness(1, 1) = 3. * bulk + elastic;
    initial_stiffness(1, 0) = initial_stiffness(0, 1) = 3. * bulk - elastic;
    initial_stiffness(2, 2) = elastic;
    initial_stiffness *= 3. * bulk / (9. * bulk - elastic);

    trial_stiffness = current_stiffness = initial_stiffness;

    initialize_history(1);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> DuncanSelig::get_copy() { return make_unique<DuncanSelig>(*this); }

int DuncanSelig::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& max_dev_stress = trial_history(0);

    // assuming elastic loading/unloading
    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    // if elastic response exceeds the maximum deviatoric stress, then the material is considered to be in a plastic state
    if(dev(trial_stress) <= max_dev_stress) return SUANPAN_SUCCESS;

    vec3 net_stress = current_stress, net_strain = incre_strain;

    if(dev(current_stress) < max_dev_stress) {
        // project onto the yield surface first
        auto elastic_portion = .5;
        if(SUANPAN_SUCCESS != project_to_surface(elastic_portion)) return SUANPAN_FAIL;

        net_stress = current_stress + initial_stiffness * incre_strain * elastic_portion;
        net_strain = (1. - elastic_portion) * incre_strain;
    }

    auto ref_error = 0.;
    vec3 incre;
    mat33 jacobian;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto [elastic, bulk, deds, dkds] = compute_moduli();

        const auto factor_a = 3. * bulk * (3. * bulk + elastic);
        const auto factor_b = 3. * bulk * (3. * bulk - elastic);

        mat33 right(fill::zeros);
        right(0, 0) = right(1, 1) = factor_a;
        right(0, 1) = right(1, 0) = factor_b;
        right(2, 2) = 3. * bulk * elastic;

        const vec3 residual = (9. * bulk - elastic) * (trial_stress - net_stress) - right * net_strain;

        const rowvec3 dfads = (18. * bulk + 3. * elastic) * dkds + 3. * bulk * deds;
        const rowvec3 dfbds = (18. * bulk - 3. * elastic) * dkds - 3. * bulk * deds;

        jacobian = (9. * bulk - elastic) * eye(3, 3) + (trial_stress - net_stress) * (9. * dkds - deds);
        jacobian.row(0) -= net_strain(0) * dfads + net_strain(1) * dfbds;
        jacobian.row(1) -= net_strain(0) * dfbds + net_strain(1) * dfads;
        jacobian.row(2) -= net_strain(2) * 3. * (bulk * deds + elastic * dkds);

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error / ref_error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) {
            if(!solve(trial_stiffness, jacobian, right)) return SUANPAN_FAIL;

            max_dev_stress = dev(trial_stress);

            return SUANPAN_SUCCESS;
        }

        trial_stress -= incre;
    }
}

int DuncanSelig::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    trial_history = current_history = initial_history;
    return SUANPAN_SUCCESS;
}

int DuncanSelig::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    current_history = trial_history;
    return SUANPAN_SUCCESS;
}

int DuncanSelig::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    trial_history = current_history;
    return SUANPAN_SUCCESS;
}

void DuncanSelig::print() {
    suanpan_info("The Duncan-Selig soil model.\n");
    suanpan_info("Strain:", current_strain);
    suanpan_info("Stress:", current_stress);
}
