/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
// ReSharper disable StringLiteralTypo
#include "DafaliasManzari.h"
#include <Toolbox/tensorToolbox.h>

const span DafaliasManzari::sb(1, 6);
const span DafaliasManzari::sk(2, 7);
const span DafaliasManzari::sl(8, 13);
const span DafaliasManzari::sm(14, 19);
const mat66 DafaliasManzari::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

DafaliasManzari::DafaliasManzari(const unsigned T, const double G0, const double NU, const double AC, const double LC, const double E0, const double XI, const double M, const double H0, const double H1, const double CH, const double NB, const double A, const double ND, const double ZM, const double CZ, const double PC, const double GR, const double R)
    : DataDafaliasManzari{fabs(G0), fabs(NU), fabs(AC), fabs(LC), fabs(E0), fabs(XI), fabs(M), fabs(H0), fabs(H1), fabs(CH), fabs(NB), A, fabs(ND), fabs(ZM), fabs(CZ), -fabs(PC), fabs(GR)}
    , Material3D(T, R) { access::rw(tolerance) = 1E-12; }

int DafaliasManzari::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(gi * (2. + 2. * poissons_ratio), poissons_ratio);

    initialize_history(18);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> DafaliasManzari::get_copy() { return make_unique<DafaliasManzari>(*this); }

double DafaliasManzari::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return gi * (2. + 2. * poissons_ratio);
    if(ParameterType::SHEARMODULUS == P || ParameterType::G == P) return gi;
    if(ParameterType::BULKMODULUS == P) return gi * (2. + 2. * poissons_ratio) / (3. - 6. * poissons_ratio);
    if(ParameterType::POISSONSRATIO == P) return poissons_ratio;
    return 0.;
}

int DafaliasManzari::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= tolerance) return SUANPAN_SUCCESS;

    const auto current_p = tensor::mean3(current_stress);
    const auto current_s = tensor::dev(current_stress);
    const auto incre_ev = tensor::trace3(incre_strain);
    const vec incre_ed = unit_dev_tensor * incre_strain;

    // assume no plasticity
    // compute trial stress

    auto p = current_p + pr * gi * incre_ev;
    vec s = current_s + 2. * gi * incre_ed;

    const auto void_ratio = e0 + (1. + e0) * tensor::trace3(trial_strain);
    const auto v_term_a = pow(2.97 - void_ratio, 2.) / (1. + void_ratio);
    const auto v_term_b = (void_ratio * (void_ratio + 2.) - 14.7609) * pow(1. + void_ratio, -2.) * (1. + e0);

    double g, pgpe, pgpp;

    vec residual(7, fill::none), incre;
    mat jacobian(7, 7, fill::eye);

    auto counter = 0u;
    auto ref_error = 1.;

    while(true) {
        if(max_iteration == ++counter) return SUANPAN_FAIL;

        const auto sqrt_term = shear_modulus * sqrt(std::max(datum::eps, pc * p));

        g = sqrt_term * v_term_a;

        if(g > gi) {
            pgpe = sqrt_term * v_term_b;
            pgpp = .5 * g / p;
        }
        else {
            g = gi;
            pgpe = pgpp = 0.;
        }

        residual(sa) = p - current_p - pr * g * incre_ev;
        residual(sb) = s - current_s - 2. * g * incre_ed;

        jacobian(sa, sa) = 1. - pr * incre_ev * pgpp;
        jacobian(sb, sa) = -2. * pgpp * incre_ed;

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        auto error = norm(residual);
        if(1 == counter) ref_error = std::max(1., error);
        suanpan_debug("DafaliasManzari local elastic iteration error: %.5E.\n", error /= ref_error);
        if(error <= tolerance) break;

        p -= incre(sa);
        s -= incre(sb);
    }

    // check if yield

    const vec current_alpha(&current_history(0), 6);

    vec eta = s + p * current_alpha;
    auto norm_eta = tensor::stress::norm(eta);

    if(norm_eta + m * p < 0.) {
        trial_stress = s + p * tensor::unit_tensor2;

        mat left(7, 6, fill::none), right;

        left.row(sa) = pr * (incre_ev * pgpe + g) * tensor::unit_tensor2.t();
        left.rows(sb) = 2. * pgpe * incre_ed * tensor::unit_tensor2.t() + 2. * g * unit_dev_tensor;

        if(!solve(right, jacobian, left)) return SUANPAN_FAIL;

        trial_stiffness = right.rows(sb);
        trial_stiffness.row(0) += right.row(sa);
        trial_stiffness.row(1) += right.row(sa);
        trial_stiffness.row(2) += right.row(sa);

        return SUANPAN_SUCCESS;
    }

    // yield function violated

    const vec current_z(&current_history(6), 6);

    trial_history = current_history;
    vec alpha(&trial_history(0), 6, false, true);
    vec z(&trial_history(6), 6, false, true);
    vec ini_alpha(&trial_history(12), 6, false, true);

    residual.set_size(20);
    jacobian.set_size(20, 20);
    jacobian(si, si) = 0.;
    jacobian(si, sm).zeros();
    jacobian(sk, sm).zeros();
    jacobian(sl, sm).zeros();

    counter = 0u;

    vec n, zz, aabmn;
    auto gamma = 0.;
    double pabpe, d, pdpe, h, phpe;
    auto update_ini_alpha = false;

    while(true) {
        if(max_iteration == ++counter) return SUANPAN_FAIL;

        // shear modulus

        auto tmp_term = shear_modulus * sqrt(std::max(datum::eps, pc * p));
        g = tmp_term * v_term_a;

        if(g > gi) {
            pgpe = tmp_term * v_term_b;
            pgpp = .5 * g / p;
        }
        else {
            g = gi;
            pgpe = pgpp = 0.;
        }

        // state parameter

        tmp_term = lc * pow(std::max(datum::eps, p / pc), xi);
        const auto psi = void_ratio - e0 + tmp_term;
        const auto ppsipp = xi * tmp_term / p;

        // surface

        const auto ad = ac * exp(nd * psi);
        const auto ab = ac * exp(-nb * psi);
        const auto adm = ad - m;
        const auto abm = ab - m;

        auto padpe = nd * ad;
        pabpe = -nb * ab;
        const auto padpp = padpe * ppsipp;
        const auto pabpp = pabpe * ppsipp;
        padpe *= 1. + e0;
        pabpe *= 1. + e0;

        // yield function

        eta = s + p * alpha;
        norm_eta = tensor::stress::norm(eta);

        n = eta / norm_eta;
        const vec unit_n = n % tensor::stress::norm_weight;
        const vec unit_alpha = alpha % tensor::stress::norm_weight;
        const auto alpha_n = dot(n, unit_alpha);
        aabmn = alpha - abm * n;

        const vec np = (alpha - alpha_n * n) / norm_eta;
        const mat ns = (eye(6, 6) - n * unit_n.t()) / norm_eta;

        // dilatancy

        const vec unit_z = z % tensor::stress::norm_weight;
        const auto zn = dot(n, unit_z);

        d = a * (adm - alpha_n);

        double pdpp;
        rowvec pdps, pdpa, pdpz;
        if(zn > 0.) {
            const auto term_a = a * (1. + zn);

            pdpe = term_a * padpe;
            pdpp = term_a * padpp + dot(d * unit_z - term_a * unit_alpha, np);
            pdps = (d * unit_z - term_a * unit_alpha).t() * ns;
            pdpa = p * pdps - term_a * unit_n.t();
            pdpz = d * unit_n.t();

            d *= 1. + zn;
        }
        else {
            pdpe = a * padpe;
            pdpp = a * (padpp - dot(unit_alpha, np));
            pdps = -a * unit_alpha.t() * ns;
            pdpa = p * pdps - a * unit_n.t();
            pdpz.zeros(6);
        }

        // hardening

        tmp_term = shear_modulus * h0 * sqrt(std::max(datum::eps, pc / p));
        const auto b0 = tmp_term * (1. - ch * void_ratio);
        const auto pb0pe = -ch * tmp_term * (1. + e0);
        const auto pb0pp = -.5 * b0 / p;

        update_ini_alpha = false;
        vec diff_alpha = (ini_alpha - alpha) % tensor::stress::norm_weight;
        tmp_term = exp(h1 * dot(diff_alpha, n));

        if(tmp_term > 1.) {
            update_ini_alpha = true;
            diff_alpha = (current_alpha - alpha) % tensor::stress::norm_weight;
            tmp_term = exp(h1 * dot(diff_alpha, n));
        }

        h = tmp_term * b0;

        phpe = tmp_term * pb0pe;
        const auto phpp = tmp_term * pb0pp;
        const rowvec phps = h * h1 * diff_alpha.t() * ns;
        const rowvec phpa = p * phps - h * h1 * unit_n.t();

        // local iteration

        residual(si) = norm_eta + m * p;
        residual(sj) = p - current_p + pr * g * (gamma * d - incre_ev);
        residual(sk) = s - current_s + 2. * g * (gamma * n - incre_ed);
        residual(sl) = alpha - current_alpha + gamma * h * aabmn;
        residual(sm) = z - current_z;

        jacobian(si, sj) = alpha_n + m;
        jacobian(si, sk) = unit_n.t();
        jacobian(si, sl) = p * jacobian(si, sk);

        const auto gk = gamma * pr * g;
        jacobian(sj, si) = d * pr * g;
        jacobian(sj, sj) = 1. + pr * pgpp * (gamma * d - incre_ev) + gk * pdpp;
        jacobian(sj, sk) = gk * pdps;
        jacobian(sj, sl) = gk * pdpa;
        jacobian(sj, sm) = gk * pdpz;

        jacobian(sk, si) = 2. * g * n;
        jacobian(sk, sj) = 2. * pgpp * (gamma * n - incre_ed) + 2. * g * gamma * np;
        jacobian(sk, sk) = 2. * g * gamma * ns;
        jacobian(sk, sl) = p * jacobian(sk, sk);
        jacobian(sk, sk) += eye(6, 6);

        jacobian(sl, si) = h * aabmn;
        jacobian(sl, sj) = gamma * phpp * aabmn - gamma * h * (pabpp * n + abm * np);
        jacobian(sl, sk) = gamma * aabmn * phps - gamma * h * abm * ns;
        jacobian(sl, sl) = (1. + gamma * h) * eye(6, 6) + gamma * aabmn * phpa - gamma * h * abm * p * ns;

        jacobian(sm, sm) = eye(6, 6);

        if(d > 0.) {
            const auto factor_a = cz * gamma;
            const auto factor_b = factor_a * d;
            const auto factor_c = factor_b * zm;

            zz = z - zm * n;

            residual(sm) += factor_b * zz;

            jacobian(sm, si) = cz * d * zz;
            jacobian(sm, sj) = factor_a * pdpp * zz - factor_c * np;
            jacobian(sm, sk) = factor_a * zz * pdps - factor_c * ns;
            jacobian(sm, sl) = factor_a * zz * pdpa - factor_c * p * ns;
            jacobian(sm, sm) += factor_b * eye(6, 6) + factor_a * zz * pdpz;
        }
        else {
            jacobian(sm, si).zeros();
            jacobian(sm, sj).zeros();
            jacobian(sm, sk).zeros();
            jacobian(sm, sl).zeros();
        }

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        auto error = norm(residual);
        if(1 == counter) ref_error = std::max(1., error);
        suanpan_debug("DafaliasManzari local plastic iteration error: %.5E.\n", error /= ref_error);
        if(error <= tolerance) break;

        gamma -= incre(si);
        p -= incre(sj);
        s -= incre(sk);
        alpha -= incre(sl);
        z -= incre(sm);
    }

    trial_stress = s + p * tensor::unit_tensor2;

    mat::fixed<20, 6> left(fill::none), right;

    left.row(si).zeros();
    left.row(sj) = pr * (pgpe * (incre_ev - gamma * d) + g - g * gamma * pdpe) * tensor::unit_tensor2.t();
    left.rows(sk) = 2. * g * unit_dev_tensor + 2. * pgpe * (incre_ed - gamma * n) * tensor::unit_tensor2.t();
    left.rows(sl) = (gamma * h * n * pabpe - gamma * aabmn * phpe) * tensor::unit_tensor2.t();

    if(d > 0.) left.rows(sm) = -cz * gamma * pdpe * zz * tensor::unit_tensor2.t();
    else left.rows(sm).zeros();

    if(!solve(right, jacobian, left)) return SUANPAN_FAIL;

    trial_stiffness = right.rows(sk);
    trial_stiffness.row(0) += right.row(sj);
    trial_stiffness.row(1) += right.row(sj);
    trial_stiffness.row(2) += right.row(sj);

    if(update_ini_alpha) ini_alpha = current_alpha;

    return SUANPAN_SUCCESS;
}

int DafaliasManzari::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int DafaliasManzari::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int DafaliasManzari::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void DafaliasManzari::print() { suanpan_info("A Dafalias--Manzari sand model. doi: 10.1061/(ASCE)0733-9399(2004)130:6(622)\n"); }
