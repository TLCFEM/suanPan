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

// ReSharper disable IdentifierTypo
#include "DuncanSelig.h"

#include <Domain/DomainBase.h>

double DuncanSelig::dev(const vec& t_stress) { return std::sqrt(std::pow(t_stress(0) - t_stress(1), 2) + std::pow(2. * t_stress(2), 2)); }

rowvec3 DuncanSelig::der_dev(const vec& t_stress) { return {t_stress(0) - t_stress(1), t_stress(1) - t_stress(0), 4. * t_stress(2)}; }

mat DuncanSelig::compute_stiffness(const double elastic, const double bulk) {
    mat stiffness(3, 3, fill::zeros);

    stiffness(0, 0) = stiffness(1, 1) = 3. * bulk + elastic;
    stiffness(1, 0) = stiffness(0, 1) = 3. * bulk - elastic;
    stiffness(2, 2) = elastic;
    stiffness *= 3. * bulk / (9. * bulk - elastic);

    return stiffness;
}

std::tuple<double, double> DuncanSelig::compute_elastic(const double s3) const {
    double elastic, deds3 = 0.;
    if(s3 < -min_ratio * p_atm) elastic = ref_elastic * std::pow(.01, n);
    else if(s3 < min_ratio * p_atm) elastic = ref_elastic * std::pow(min_ratio, n);
    else deds3 = n * (elastic = ref_elastic * std::pow(s3 / p_atm, n)) / s3;

    return {elastic, deds3};
}

std::tuple<double, double> DuncanSelig::compute_bulk(const double s3) const {
    double bulk, dkds3 = 0.;
    if(s3 < -min_ratio * p_atm) bulk = ref_bulk * std::pow(.01, m);
    else if(s3 < min_ratio * p_atm) bulk = ref_bulk * std::pow(min_ratio, m);
    else dkds3 = m * (bulk = ref_bulk * std::pow(s3 / p_atm, m)) / s3;

    return {bulk, dkds3};
}

DuncanSelig::ds_moduli DuncanSelig::compute_elastic_moduli() {
    // principal stresses

    const auto radius = .5 * dev(trial_stress);
    rowvec3 drds(fill::zeros);
    if(radius > datum::eps) drds = der_dev(trial_stress) / radius * .25;

    const auto center = -.5 * (trial_stress(0) + trial_stress(1));
    const rowvec3 dcds{-.5, -.5, 0.};

    const auto s3 = center - radius;

    const rowvec3 ds3ds = dcds - drds;

    // for elastic modulus

    const auto [elastic, deds3] = compute_elastic(s3);

    const rowvec3 deds = deds3 * ds3ds;

    // for bulk modulus

    auto [bulk, dkds3] = compute_bulk(s3);

    rowvec3 dkds = dkds3 * ds3ds;

    if(3. * bulk < elastic) {
        bulk = elastic / 3.;
        dkds = deds / 3.;
    }

    return {elastic, bulk, deds, dkds};
}

DuncanSelig::ds_moduli DuncanSelig::compute_plastic_moduli() {
    // principal stresses

    const auto radius = .5 * dev(trial_stress);
    rowvec3 drds(fill::zeros);
    if(radius > datum::eps) drds = der_dev(trial_stress) / radius * .25;

    const auto center = -.5 * (trial_stress(0) + trial_stress(1));
    const rowvec3 dcds{-.5, -.5, 0.};

    const auto s1 = center + radius, s3 = center - radius;

    const rowvec3 ds1ds = dcds + drds, ds3ds = dcds - drds;

    // for elastic modulus

    double phi, dphids3 = 0.;
    if(s3 < p_atm) phi = ini_phi;
    else if(phi = ini_phi - ten_fold_phi_diff * log10(s3 / p_atm); phi < 0.) phi = 0.;
    else dphids3 = -ten_fold_phi_diff / (s3 * log(10));

    const auto denom = 1. - std::sin(phi);
    auto max_dev_stress = 2. / r_f * (cohesion * std::cos(phi) + s3 * std::sin(phi)) / denom;
    auto dmdsds3 = 0.;
    if(max_dev_stress > min_ratio * p_atm) {
        const auto pmdspphi = 2. / r_f * (s3 * std::cos(phi) / denom + cohesion) / denom;
        const auto pmdsps3 = 2. / r_f * std::sin(phi) / denom;
        dmdsds3 = pmdspphi * dphids3 + pmdsps3;
    }
    else max_dev_stress = min_ratio * p_atm;

    const auto dev_stress = s1 - s3;
    const auto pdsps1 = 1.;
    const auto pdsps3 = -1.;

    const auto [ini_elastic, deids3] = compute_elastic(s3);

    const auto pepei = std::pow(1. - dev_stress / max_dev_stress, 2.);
    const auto elastic = ini_elastic * pepei;
    const auto pepds = -2. * ini_elastic * (1. - dev_stress / max_dev_stress) / max_dev_stress;
    const auto pepmds = -pepds * dev_stress / max_dev_stress;

    const auto peps1 = pepds * pdsps1;
    const auto peps3 = pepei * deids3 + pepds * pdsps3 + pepmds * dmdsds3;

    const rowvec3 deds = peps1 * ds1ds + peps3 * ds3ds;

    // for bulk modulus

    auto [bulk, dkds3] = compute_bulk(s3);

    rowvec3 dkds = dkds3 * ds3ds;

    if(3. * bulk < elastic) {
        bulk = elastic / 3.;
        dkds = deds / 3.;
    }

    return {elastic, bulk, deds, dkds};
}

int DuncanSelig::project_onto_surface(double& multiplier) {
    const auto max_dev_stress = trial_history(0);

    trial_stress = current_stress + current_stiffness * incre_strain;
    multiplier = 1.;

    auto ref_error = 0.;
    vec4 incre, residual;

    static const uvec sa{0, 1, 2}, sb{3};

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Local elastic iteration cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto [elastic, bulk, deds, dkds] = compute_elastic_moduli();

        const auto factor_a = 3. * bulk * (3. * bulk + elastic);
        const auto factor_b = 3. * bulk * (3. * bulk - elastic);

        mat33 right(fill::zeros);
        right(0, 0) = right(1, 1) = factor_a;
        right(0, 1) = right(1, 0) = factor_b;
        right(2, 2) = 3. * bulk * elastic;

        const auto trial_dev_stress = dev(trial_stress);
        const vec3 t_stress = trial_stress - current_stress;
        const vec3 t_strain = incre_strain * multiplier;

        residual.head(3) = (9. * bulk - elastic) * t_stress - right * t_strain;
        residual(3) = trial_dev_stress - max_dev_stress;

        const rowvec3 factor_c = 3. * (elastic * dkds + bulk * deds);
        const rowvec3 dfads = 18. * bulk * dkds + factor_c;
        const rowvec3 dfbds = 18. * bulk * dkds - factor_c;

        mat44 jacobian(fill::zeros);
        jacobian(sa, sa) = (9. * bulk - elastic) * eye(3, 3) + t_stress * (9. * dkds - deds);
        jacobian.row(0).head(3) -= t_strain(0) * dfads + t_strain(1) * dfbds;
        jacobian.row(1).head(3) -= t_strain(0) * dfbds + t_strain(1) * dfads;
        jacobian.row(2).head(3) -= t_strain(2) * factor_c;
        jacobian(sa, sb) = -right * incre_strain;
        if(trial_dev_stress > datum::eps) jacobian(sb, sa) = der_dev(trial_stress) / trial_dev_stress;

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local elastic iteration error: {:.5E}.\n", error / ref_error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) return SUANPAN_SUCCESS;

        trial_stress -= incre.head(3);
        multiplier -= incre(3);
    }
}

int DuncanSelig::local_update(const vec& ref_stress, const vec& ref_strain, const bool two_stage) {
    const auto update_moduli = [&] {
        const auto max_dev_stress = trial_history(0);
        return two_stage || dev(trial_stress) > max_dev_stress ? compute_plastic_moduli() : compute_elastic_moduli();
    };

    auto ref_error = 0.;
    vec3 incre;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            if(two_stage)
                suanpan_error("Local iteration cannot converge within {} iterations.\n", max_iteration);
            else
                suanpan_debug("Local iteration cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto [elastic, bulk, deds, dkds] = update_moduli();

        const auto factor_a = 3. * bulk * (3. * bulk + elastic);
        const auto factor_b = 3. * bulk * (3. * bulk - elastic);

        mat33 right(fill::zeros);
        right(0, 0) = right(1, 1) = factor_a;
        right(0, 1) = right(1, 0) = factor_b;
        right(2, 2) = 3. * bulk * elastic;

        const vec3 t_stress = trial_stress - ref_stress;
        const vec3 residual = (9. * bulk - elastic) * t_stress - right * ref_strain;

        const rowvec3 factor_c = 3. * (elastic * dkds + bulk * deds);
        const rowvec3 dfads = 18. * bulk * dkds + factor_c;
        const rowvec3 dfbds = 18. * bulk * dkds - factor_c;

        mat33 jacobian = (9. * bulk - elastic) * eye(3, 3) + t_stress * (9. * dkds - deds);
        jacobian.row(0) -= ref_strain(0) * dfads + ref_strain(1) * dfbds;
        jacobian.row(1) -= ref_strain(0) * dfbds + ref_strain(1) * dfads;
        jacobian.row(2) -= ref_strain(2) * factor_c;

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) {
            if(!solve(trial_stiffness, jacobian, right)) return SUANPAN_FAIL;

            return SUANPAN_SUCCESS;
        }

        trial_stress -= incre;
    }
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
    const auto [elastic, bulk, deds, dkds] = compute_elastic_moduli();

    trial_stiffness = current_stiffness = initial_stiffness = compute_stiffness(elastic, bulk);

    initialize_history(1);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> DuncanSelig::get_copy() { return std::make_unique<DuncanSelig>(*this); }

int DuncanSelig::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;

    trial_stress = current_stress + current_stiffness * incre_strain;

    const auto update_dev_stress = [&] {
        auto& max_dev_stress = trial_history(0);
        if(const auto trail_dev_stress = dev(trial_stress); trail_dev_stress > max_dev_stress) max_dev_stress = trail_dev_stress;
        return SUANPAN_SUCCESS;
    };

    // first try a whole step size local iteration
    if(SUANPAN_SUCCESS == local_update(current_stress, incre_strain, false)) return update_dev_stress();

    // alternatively, directly use the plastic branch
    // if(SUANPAN_SUCCESS == local_update(current_stress, incre_strain, true)) return update_dev_stress();

    // if that fails, mostly likely due to discontinuity of the gradient
    // assume proportional loading, project the stress onto the yield surface using elastic moduli
    auto multiplier = 1.;
    if(SUANPAN_SUCCESS != project_onto_surface(multiplier)) return SUANPAN_FAIL;

    // then using plastic moduli to compute the new plastic state
    // !!! the tangent operator is not algorithmically consistent in this case
    if(SUANPAN_SUCCESS != local_update(vec(trial_stress), (1. - multiplier) * incre_strain, true)) return SUANPAN_FAIL;

    return update_dev_stress();
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
    suanpan_info("The Duncan--Selig soil model.\n");
    suanpan_info("Strain:", current_strain);
    suanpan_info("Stress:", current_stress);
}
