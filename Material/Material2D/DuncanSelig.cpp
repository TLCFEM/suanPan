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

std::tuple<double, double> DuncanSelig::compute_modulus(const double ps, const double ref_modulus, const double order) const {
    const auto fraction = ps / p_atm;
    if(fraction <= min_ratio) return {ref_modulus * std::pow(min_ratio, order), 0.};

    const auto modulus = ref_modulus * std::pow(fraction, order);
    return {modulus, order * modulus / ps};
}

std::tuple<double, double> DuncanSelig::compute_elastic(const double ps) const { return compute_modulus(ps, ref_elastic, n); }

std::tuple<double, double> DuncanSelig::compute_bulk(const double ps) const { return compute_modulus(ps, ref_bulk, m); }

DuncanSelig::ds_tuple DuncanSelig::compute_principal_stress(const vec& t_stress) const {
    constexpr auto center_scale{-.5}; // this is correct as later the algorithm uses s3
    const auto center = center_scale * (t_stress(0) + t_stress(1));
    const rowvec3 dcds{center_scale, center_scale, 0.};

    const auto radius = .5 * dev(t_stress);
    rowvec3 drds(fill::zeros);
    if(radius > tolerance * std::max(1., std::fabs(center))) drds = .25 / radius * rowvec3{t_stress(0) - t_stress(1), t_stress(1) - t_stress(0), 4. * t_stress(2)};

    return {center + radius, center - radius, dcds + drds, dcds - drds};
}

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

mat DuncanSelig::compute_stiffness(const vec& t_stress) {
    const auto [elastic, bulk, deds, dkds] = compute_elastic_moduli(t_stress);

    mat stiffness(3, 3, fill::zeros);

    stiffness(0, 0) = stiffness(1, 1) = 3. * bulk + elastic;
    stiffness(1, 0) = stiffness(0, 1) = 3. * bulk - elastic;
    stiffness(2, 2) = elastic;
    stiffness *= 3. * bulk / (9. * bulk - elastic);

    return stiffness;
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
    trial_stiffness = current_stiffness = initial_stiffness = compute_stiffness(vec{3, fill::zeros});

    initialize_history(1);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> DuncanSelig::unique_copy() { return std::make_unique<DuncanSelig>(*this); }

int DuncanSelig::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    trial_stress = current_stress + compute_stiffness(current_stress) * incre_strain;

    // first assume elastic loading/unloading
    if(SUANPAN_SUCCESS != local_update<false>()) return SUANPAN_FAIL;

    // if the elastic solution violates the plastic condition, redo the local iteration with plastic moduli
    if(dev(trial_stress) > trial_history(0) && SUANPAN_SUCCESS != local_update<true>()) return SUANPAN_FAIL;

    trial_history(0) = std::max(trial_history(0), dev(trial_stress));

    return SUANPAN_SUCCESS;
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
