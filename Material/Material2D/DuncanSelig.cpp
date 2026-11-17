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

// ReSharper disable IdentifierTypo
#include "DuncanSelig.h"

#include <Domain/DomainBase.h>

double DuncanSelig::dev(const vec& t_stress) { /* omits the 0.5 that is moved to other places */ return std::hypot(t_stress(0) - t_stress(1), 2. * t_stress(2)); }

mat DuncanSelig::compute_stiffness(const double elastic, const double bulk) {
    mat stiffness(3, 3, fill::zeros);

    stiffness(0, 0) = stiffness(1, 1) = 3. * bulk + elastic;
    stiffness(1, 0) = stiffness(0, 1) = 3. * bulk - elastic;
    stiffness(2, 2) = elastic;
    stiffness *= 3. * bulk / (9. * bulk - elastic);

    return stiffness;
}

DuncanSelig::ds_tuple DuncanSelig::compute_principal_stress(const vec& t_stress) {
    const auto radius = .5 * dev(t_stress);

    static constexpr auto center_scale{-.5}; // this is correct as later the algorithm uses s3
    const auto center = center_scale * (t_stress(0) + t_stress(1));
    const rowvec3 dcds{center_scale, center_scale, 0.};

    rowvec3 drds(fill::zeros);
    // the gradient of the radius is not unique at a hydrostatic stress state, in which the
    // deviatoric stress is dominated by floating-point round-off, use the isotropic gradient
    if(radius > 1E-10 * std::max(1., std::fabs(center))) drds = .25 / radius * rowvec3{t_stress(0) - t_stress(1), t_stress(1) - t_stress(0), 4. * t_stress(2)};

    return {center + radius, center - radius, dcds + drds, dcds - drds};
}

std::tuple<double, double> DuncanSelig::compute_modulus(const double ps, const double ref_modulus, const double order) const {
    const auto fraction = ps / p_atm;
    if(fraction <= min_ratio) return {ref_modulus * std::pow(min_ratio, order), 0.};

    const auto modulus = ref_modulus * std::pow(fraction, order);
    return {modulus, order * modulus / ps};
}

std::tuple<double, double> DuncanSelig::compute_elastic(const double ps) const { return compute_modulus(ps, ref_elastic, n); }

std::tuple<double, double> DuncanSelig::compute_bulk(const double ps) const { return compute_modulus(ps, ref_bulk, m); }

DuncanSelig::ds_tuple DuncanSelig::compute_elastic_moduli(const vec& t_stress) const {
    // principal stresses

    const auto [s1, s3, ds1ds, ds3ds] = compute_principal_stress(t_stress);

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

DuncanSelig::ds_tuple DuncanSelig::compute_plastic_moduli(const vec& t_stress) const {
    // principal stresses

    const auto [s1, s3, ds1ds, ds3ds] = compute_principal_stress(t_stress);

    // for elastic modulus

    double phi, dphids3 = 0.;
    if(s3 < p_atm) phi = ini_phi;
    else if(phi = ini_phi - ten_fold_phi_diff * log10(s3 / p_atm); phi < 0.) phi = 0.;
    else dphids3 = -ten_fold_phi_diff / (s3 * log(10));

    const auto denom = 1. - std::sin(phi);
    auto max_dev_stress = 2. / r_f * (cohesion * std::cos(phi) + std::max(s3, 0.) * std::sin(phi)) / denom;
    auto dmdsds3 = 0.;
    if(max_dev_stress <= min_ratio * p_atm) max_dev_stress = min_ratio * p_atm;
    else if(s3 > 0.) {
        const auto pmdspphi = 2. / r_f * (s3 * std::cos(phi) / denom + cohesion) / denom;
        const auto pmdsps3 = 2. / r_f * std::sin(phi) / denom;
        dmdsds3 = pmdspphi * dphids3 + pmdsps3;
    }
    else dmdsds3 = 2. / r_f * cohesion / denom * dphids3;

    const auto dev_stress = s1 - s3;
    constexpr auto pdsps1 = 1.;
    constexpr auto pdsps3 = -1.;

    const auto [ini_elastic, deids3] = compute_elastic(s3);

    auto dddds{0.}, dddmds{0.};
    auto degradation = 1. - dev_stress / max_dev_stress;
    if(degradation < min_degradation) degradation = min_degradation;
    else {
        dddds = -1. / max_dev_stress;
        dddmds = dev_stress / max_dev_stress / max_dev_stress;
    }
    const auto pepei = degradation * degradation;
    const auto elastic = ini_elastic * pepei;
    const auto pepds = 2. * ini_elastic * degradation * dddds;
    const auto pepmds = 2. * ini_elastic * degradation * dddmds;

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
    const auto [elastic, bulk, deds, dkds] = compute_elastic_moduli(vec{3, fill::zeros});

    trial_stiffness = current_stiffness = initial_stiffness = compute_stiffness(elastic, bulk);

    initialize_history(1);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> DuncanSelig::unique_copy() { return std::make_unique<DuncanSelig>(*this); }

int DuncanSelig::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + current_stiffness * incre_strain;
    trial_history = current_history;

    // improve the predictor by a few fixed-point iterations so that large stress increments do not
    // start the local Newton iteration from a stress state that is far from the converged solution
    for(auto I = 0u; I < 3u; ++I) {
        const auto moduli = dev(trial_stress) > trial_history(0) ? compute_plastic_moduli(trial_stress) : compute_elastic_moduli(trial_stress);
        trial_stress = current_stress + compute_stiffness(std::get<0>(moduli), std::get<1>(moduli)) * incre_strain;
    }

    auto ref_error = 0.;
    vec3 incre;

    // the plastic/elastic branch is discontinuous at the historical maximum deviatoric stress,
    // switching branches within the local Newton iteration can result in an endless oscillation
    // between the two branches, freeze the branch determined at the predictor for this increment
    const auto plastic = dev(trial_stress) > trial_history(0);

    auto counter = 0u;
    while(true) {
        if(max_iteration < ++counter) {
            suanpan_debug("Local iteration cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto double_dev_stress = dev(trial_stress);

        const auto [elastic, bulk, deds, dkds] = plastic ? compute_plastic_moduli(trial_stress) : compute_elastic_moduli(trial_stress);

        const auto factor_a = 3. * bulk * (3. * bulk + elastic);
        const auto factor_b = 3. * bulk * (3. * bulk - elastic);

        mat33 right(fill::zeros);
        right(0, 0) = right(1, 1) = factor_a;
        right(0, 1) = right(1, 0) = factor_b;
        right(2, 2) = 3. * bulk * elastic;

        const vec3 t_stress = trial_stress - current_stress;
        const vec3 residual = (9. * bulk - elastic) * t_stress - right * incre_strain;

        const rowvec3 factor_c = 3. * (elastic * dkds + bulk * deds);
        const rowvec3 dfads = 18. * bulk * dkds + factor_c;
        const rowvec3 dfbds = 18. * bulk * dkds - factor_c;

        mat33 jacobian = (9. * bulk - elastic) * eye(3, 3) + t_stress * (9. * dkds - deds);
        jacobian.row(0) -= incre_strain(0) * dfads + incre_strain(1) * dfbds;
        jacobian.row(1) -= incre_strain(0) * dfbds + incre_strain(1) * dfads;
        jacobian.row(2) -= incre_strain(2) * factor_c;

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = suanpan::inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || suanpan::inf_norm(residual) < tolerance) && counter > 5u)) {
            if(!solve(trial_stiffness, jacobian, right)) return SUANPAN_FAIL;

            trial_history(0) = std::max(trial_history(0), double_dev_stress);
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
    suanpan_info("The Duncan--Selig soil model.\n");
    suanpan_info("Strain:", current_strain);
    suanpan_info("Stress:", current_stress);
}
