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
#include "CDPM2.h"

#include <Recorder/OutputType.h>
#include <Toolbox/ridders.hpp>
#include <Toolbox/tensor.h>
#include <Toolbox/utility.h>

const double CDPM2::sqrt_six = std::sqrt(6.);
const double CDPM2::sqrt_three_two = std::sqrt(1.5);
const mat CDPM2::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

void CDPM2::compute_plasticity(const double lode, const double s, const double p, const double kp, vec18& data) const {
    auto& f = data(0);
    auto& pfps = data(1);
    auto& pfpp = data(2);
    auto& pfpkp = data(3);
    auto& gs = data(4);
    auto& gp = data(5);
    auto& gg = data(6);
    auto& pgsps = data(7);
    auto& pgspp = data(8);
    auto& pgspkp = data(9);
    auto& pgpps = data(10);
    auto& pgppp = data(11);
    auto& pgppkp = data(12);
    auto& xh = data(13);
    auto& dxhdp = data(14);
    auto& pfpl = data(15);
    auto& r = data(16);
    auto& drdl = data(17);

    auto qh1 = 1., qh2 = 1.;
    auto dqh1dkp = 0., dqh2dkp = 0.;

    if(kp < 1.) {
        qh1 = qh0 + (1. - qh0) * kp * (kp * (kp - 3.) + 3.) - hp * kp * (kp - 1.) * (kp - 2.);
        dqh1dkp = (3. - 3. * qh0) * std::pow(kp - 1., 2.) - hp * (3. * kp * (kp - 2.) + 2.);
    }
    else {
        qh2 = 1. + hp * (kp - 1.);
        dqh2dkp = hp;
    }

    const auto ag = 3. * ftfc * qh2 + .5 * m0;
    const auto dagdkp = 3. * ftfc * dqh2dkp;
    const auto cg = qh2 * (1. + ftfc) / 3.;
    const auto dg = std::log(ag) - std::log(3. * qh2 + .5 * m0) + lndf;
    const auto dcgdkp = dqh2dkp * (1. + ftfc) / 3.;
    const auto ddgdkp = dagdkp / ag - dqh2dkp / (qh2 + m0 / 6.);
    const auto bg = cg / dg;
    const auto dbgdkp = (dcgdkp - bg * ddgdkp) / dg;

    const auto eg = (p / fc - qh2 * ftfc / 3.) / bg;
    const auto pegpkp = (-ftfc / 3. * dqh2dkp - eg * dbgdkp) / bg;
    const auto pmgpp = ag * std::exp(eg);

    const auto g3 = (s / sqrt_six + p) / fc;
    const auto g1 = (1. - qh1) * g3 * g3 + sqrt_three_two * s / fc;

    const auto square_term = ra * lode * lode;
    const auto sqrt_term = std::sqrt(rb * square_term + rc);
    const auto numerator = square_term + rb;
    const auto denominator = ra * lode + sqrt_term;
    r = numerator / denominator;
    drdl = ra / denominator * (2. * lode - r - r * rb * lode / sqrt_term);

    const auto g4 = (r * s / sqrt_six + p) / fc;
    const auto pg4pp = 1. / fc;
    const auto pg4ps = r / sqrt_six / fc;
    const auto pg4pl = s / sqrt_six / fc * drdl;

    const auto pg3pp = 1. / fc;
    const auto pg3ps = pg3pp / sqrt_six;

    const auto pg2pp = pmgpp / fc;
    const auto pg2ps = m0 * pg3ps;

    const auto pg1pp = (2. - 2. * qh1) * g3 * pg3pp;
    const auto pg1ps = (2. - 2. * qh1) * g3 * pg3ps + sqrt_three_two / fc;
    const auto pg1pkp = -dqh1dkp * g3 * g3;

    f = g1 * g1 + m0 * qh1 * qh1 * qh2 * g4 - qh1 * qh1 * qh2 * qh2;

    pfpp = 2. * g1 * pg1pp + m0 * qh1 * qh1 * qh2 * pg4pp;
    pfps = 2. * g1 * pg1ps + m0 * qh1 * qh1 * qh2 * pg4ps;
    pfpkp = 2. * g1 * pg1pkp + 2. * qh1 * qh2 * (m0 * g4 * dqh1dkp - qh1 * dqh2dkp - qh2 * dqh1dkp) + m0 * qh1 * qh1 * g4 * dqh2dkp;
    pfpl = m0 * qh1 * qh1 * qh2 * pg4pl;

    gp = 2. * g1 * pg1pp + qh1 * qh1 * pg2pp;
    gs = 2. * g1 * pg1ps + qh1 * qh1 * pg2ps;

    gg = std::sqrt(gs * gs + gp * gp / 3.);

    pgppp = (4. - 4. * qh1) * (pg1pp * g3 + g1 * pg3pp) + qh1 * qh1 * ag * std::exp(eg) / bg / fc;
    pgpps = (4. - 4. * qh1) * (pg1ps * g3 + g1 * pg3ps);
    pgppkp = 4. * g3 * (pg1pkp - qh1 * pg1pkp - dqh1dkp * g1) + (2. * dqh1dkp * ag + qh1 * (dagdkp + ag * pegpkp)) * qh1 * std::exp(eg);

    pgppp /= fc;
    pgpps /= fc;
    pgppkp /= fc;

    pgspp = (4. - 4. * qh1) * (pg1pp * g3 + g1 * pg3pp) + 6. * pg1pp;
    pgsps = (4. - 4. * qh1) * (pg1ps * g3 + g1 * pg3ps) + 6. * pg1ps;
    pgspkp = 4. * g3 * (pg1pkp - qh1 * pg1pkp - dqh1dkp * g1) + 6. * pg1pkp + 2. * m0 * qh1 * dqh1dkp;

    pgspp /= sqrt_six * fc;
    pgsps /= sqrt_six * fc;
    pgspkp /= sqrt_six * fc;

    if(const auto rh = -p / fc - 1. / 3.; rh >= 0.) {
        xh = (bh - ah) * std::exp(-rh / ch);
        dxhdp = xh / ch / fc;
        xh += ah;
    }
    else {
        xh = eh * std::exp(rh / fh);
        dxhdp = -xh / fh / fc;
        xh += dh;
    }
}

int CDPM2::compute_damage(const double gamma, const double s, const double p, const double kp, const double ac, vec18& data) {
    const auto gs = data(4);
    const auto gp = data(5);
    const auto gg = data(6);
    const auto pgsps = data(7);
    const auto pgspp = data(8);
    const auto pgspkp = data(9);
    const auto pgpps = data(10);
    const auto pgppp = data(11);
    const auto pgppkp = data(12);
    const auto r = data(16);
    const auto drdl = data(17);

    const auto& current_ee = current_history(7);
    const auto& current_et = current_history(8);
    const auto& current_ec = current_history(9);
    const auto& current_kdt1 = current_history(12);
    const auto& current_kdc1 = current_history(13);
    const auto& current_kdt2 = current_history(14);
    const auto& current_kdc2 = current_history(15);
    auto& ee = trial_history(7);
    auto& et = trial_history(8);
    auto& ec = trial_history(9);
    auto& kdt = trial_history(10);
    auto& kdc = trial_history(11);
    auto& kdt1 = trial_history(12);
    auto& kdc1 = trial_history(13);
    auto& kdt2 = trial_history(14);
    auto& kdc2 = trial_history(15);
    auto& omegat = trial_history(16);
    auto& omegac = trial_history(17);

    // ee
    const auto ptapp = .5 * e0 * m0 / fc;
    const auto ptaps = ptapp / sqrt_six * r;
    const auto ptapl = ptapp / sqrt_six * s * drdl;
    const auto ptbps = sqrt_three_two * e0 / fc;
    const auto term_a = ptaps * s + ptapp * p;
    const auto term_b = ptbps * s;
    const auto term_c = std::sqrt(term_a * term_a + term_b * term_b);
    ee = term_a + term_c;
    const auto incre_ee = ee - current_ee;
    const auto peeps = (ptaps * ee + term_b * ptbps) / term_c;
    const auto peepp = ptapp * ee / term_c;
    const auto peepl = ptapl * ee / term_c;

    // ep
    const auto ep = gamma * gg;
    const auto peppg = gg;
    const auto pepps = gamma / gg * (gs * pgsps + gp / 3. * pgpps);
    const auto peppp = gamma / gg * (gs * pgspp + gp / 3. * pgppp);
    const auto peppkp = gamma / gg * (gs * pgspkp + gp / 3. * pgppkp);

    // xs
    auto xs = 1., pxsps = 0., pxspp = 0.;
    if(p <= 0.) {
        pxspp = (sqrt_six - as * sqrt_six) / s;
        xs += pxspp * p;
        pxsps = -pxspp * p / s;
    }

    const auto at = 1. - ac;

    // kdt
    auto incre_kdt = 0., pkdtps = 0., pkdtpp = 0., pkdtpac = 0., pkdtpl = 0.;
    if((et = current_et + at * incre_ee) > kdt) {
        incre_kdt = et - kdt;
        kdt = et;
        pkdtps = at * peeps;
        pkdtpp = at * peepp;
        pkdtpac = -incre_ee;
        pkdtpl = at * peepl;
    }

    // kdt1
    auto pkdt1pg = 0., pkdt1ps = 0., pkdt1pp = 0., pkdt1pkp = 0., pkdt1pac = 0.;
    if(incre_kdt > 0. && kdt > e0) {
        const auto atxs = at / xs;
        const auto incre_kdt1 = atxs * ep;
        kdt1 = current_kdt1 + incre_kdt1;
        pkdt1pg = atxs * peppg;
        pkdt1ps = atxs * (pepps - incre_kdt1 * pxsps);
        pkdt1pp = atxs * (peppp - incre_kdt1 * pxspp);
        pkdt1pkp = atxs * peppkp;
        pkdt1pac = -ep / xs;
    }

    // kdt2
    const auto incre_kdt2 = incre_kdt / xs;
    kdt2 = current_kdt2 + incre_kdt2;
    const auto pkdt2ps = (pkdtps - incre_kdt2 * pxsps) / xs;
    const auto pkdt2pp = (pkdtpp - incre_kdt2 * pxspp) / xs;
    const auto pkdt2pac = pkdtpac / xs;
    const auto pkdt2pl = pkdtpl / xs;

    // kdc
    auto incre_kdc = 0., pkdcps = 0., pkdcpp = 0., pkdcpac = 0., pkdcpl = 0.;
    if((ec = current_ec + ac * incre_ee) > kdc) {
        incre_kdc = ec - kdc;
        kdc = ec;
        pkdcps = ac * peeps;
        pkdcpp = ac * peepp;
        pkdcpac = incre_ee;
        pkdcpl = ac * peepl;
    }

    // kdc1
    auto pkdc1pg = 0., pkdc1ps = 0., pkdc1pp = 0., pkdc1pkp = 0., pkdc1pac = 0.;
    if(incre_kdc > 0. && kdc > e0) {
        auto qh2 = 1., dqh2dkp = 0.;
        if(kp >= 1.) {
            qh2 += hp * kp - hp;
            dqh2dkp = hp;
        }

        const auto betac = sqrtdf * qh2 / s;
        const auto pbetacpkp = sqrtdf / s * dqh2dkp;
        const auto pbetacps = -betac / s;

        pkdc1pac = ep * betac / xs;
        const auto incre_kdc1 = pkdc1pac * ac;
        kdc1 = current_kdc1 + incre_kdc1;
        pkdc1pg = peppg * ac * betac / xs;
        pkdc1ps = ac / xs * (pepps * betac + ep * pbetacps - ep * betac / xs * pxsps);
        pkdc1pp = ac / xs * betac * (peppp - ep / xs * pxspp);
        pkdc1pkp = ac / xs * (peppkp * betac + ep * pbetacpkp);
    }

    // kdc2
    const auto incre_kdc2 = incre_kdc / xs;
    kdc2 = current_kdc2 + incre_kdc2;
    const auto pkdc2ps = (pkdcps - incre_kdc2 * pxsps) / xs;
    const auto pkdc2pp = (pkdcpp - incre_kdc2 * pxspp) / xs;
    const auto pkdc2pac = pkdcpac / xs;
    const auto pkdc2pl = pkdcpl / xs;

    vec datad(3);

    if(SUANPAN_SUCCESS != compute_damage_factor(kdt, kdt1, kdt2, eft, omegat, datad)) return SUANPAN_FAIL;
    const auto& potpkdt = datad(0);
    const auto& potpkdt1 = datad(1);
    const auto& potpkdt2 = datad(2);

    auto& potpg = data(0);
    auto& potps = data(1);
    auto& potpq = data(2);
    auto& potpkp = data(3);
    auto& potpac = data(8);
    auto& potpl = data(9);

    potpg = potpkdt1 * pkdt1pg;
    potps = potpkdt * pkdtps + potpkdt1 * pkdt1ps + potpkdt2 * pkdt2ps;
    potpq = potpkdt * pkdtpp + potpkdt1 * pkdt1pp + potpkdt2 * pkdt2pp;
    potpkp = potpkdt1 * pkdt1pkp;
    potpac = potpkdt * pkdtpac + potpkdt1 * pkdt1pac + potpkdt2 * pkdt2pac;
    potpl = potpkdt * pkdtpl + potpkdt2 * pkdt2pl;

    if(SUANPAN_SUCCESS != compute_damage_factor(kdc, kdc1, kdc2, efc, omegac, datad)) return SUANPAN_FAIL;
    const auto& pocpkdc = datad(0);
    const auto& pocpkdc1 = datad(1);
    const auto& pocpkdc2 = datad(2);

    auto& pocpg = data(4);
    auto& pocps = data(5);
    auto& pocpq = data(6);
    auto& pocpkp = data(7);
    auto& pocpac = data(10);
    auto& pocpl = data(11);

    pocpg = pocpkdc1 * pkdc1pg;
    pocps = pocpkdc * pkdcps + pocpkdc1 * pkdc1ps + pocpkdc2 * pkdc2ps;
    pocpq = pocpkdc * pkdcpp + pocpkdc1 * pkdc1pp + pocpkdc2 * pkdc2pp;
    pocpkp = pocpkdc1 * pkdc1pkp;
    pocpac = pocpkdc * pkdcpac + pocpkdc1 * pkdc1pac + pocpkdc2 * pkdc2pac;
    pocpl = pocpkdc * pkdcpl + pocpkdc2 * pkdc2pl;

    return SUANPAN_SUCCESS;
}

int CDPM2::compute_damage_factor(const double kd, const double kd1, const double kd2, const double ef, double& omega, vec& data) const {
    auto& popkd = data(0);
    auto& popkd1 = data(1);
    auto& popkd2 = data(2);

    if(kd < e0) {
        popkd = 0.;
        popkd1 = 0.;
        popkd2 = 0.;
        return SUANPAN_SUCCESS;
    }

    omega = 1.; // initial guess with the maximum value to avoid convergence to the negative solution

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) return SUANPAN_FAIL;

        const auto term_a = ft * std::exp(-(kd1 + omega * kd2) / ef);
        const auto term_b = (1. - omega) * elastic_modulus;
        const auto residual = term_a - term_b * kd;
        const auto jacobian = elastic_modulus * kd - kd2 / ef * term_a;
        const auto incre = residual / jacobian;

        const auto error = std::fabs(incre);
        suanpan_debug("Local damage iteration error: {:.5E}.\n", error);

        if(error < tolerance || (std::fabs(residual) < tolerance && counter > 5u)) {
            popkd = term_b / jacobian;
            popkd1 = term_a / ef / jacobian;
            popkd2 = popkd1 * omega;
            return SUANPAN_SUCCESS;
        }

        omega -= incre;
    }
}

CDPM2::CDPM2(const unsigned T, const double E, const double V, const double FT, const double FC, const double QH0, const double HP, const double DF, const double AH, const double BH, const double CH, const double DH, const double AS, const double EFT, const double EFC, const DamageType DT, const double R)
    : DataCDPM2{std::fabs(E), std::fabs(V), std::fabs(FT), std::fabs(FC), std::fabs(QH0), std::max(HP, static_cast<double>(std::numeric_limits<float>::epsilon())), DF, AH, BH, CH, DH, AS, std::fabs(EFT), std::fabs(EFC)}
    , Material3D(T, R)
    , damage_type(DT) {}

int CDPM2::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(18);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> CDPM2::get_copy() { return std::make_unique<CDPM2>(*this); }

double CDPM2::get_parameter(const ParameterType P) const { return material_property(elastic_modulus, poissons_ratio)(P); }

int CDPM2::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    const auto& current_kp = current_history(0);
    auto& kp = trial_history(0);
    vec plastic_strain(&trial_history(1), 6, false, true);

    trial_stress = (trial_stiffness = initial_stiffness) * (trial_strain - plastic_strain);

    //
    // plasticity part
    //

    const auto dev_stress = tensor::dev(trial_stress);
    const auto hydro_stress = tensor::mean3(trial_stress);
    const auto trial_s = tensor::stress::norm(dev_stress);
    const auto trial_p = hydro_stress;
    const vec n = dev_stress / trial_s;

    static constexpr double low_limit = -.95;
    static const double low_slope = (2. * std::cos(std::acos(low_limit) / 3.) - 1.) / (low_limit + 1.);
    static constexpr double high_limit = .95;
    static const double high_slope = (2. * std::cos(std::acos(high_limit) / 3.) - 2.) / (high_limit - 1.);

    double lode, dlode;
    if(const auto lode_a = tensor::stress::lode(dev_stress); lode_a < low_limit) {
        // close to left boundary
        // use linear approximation
        lode = 1. + low_slope * (1. + lode_a);
        dlode = low_slope;
    }
    else if(lode_a > high_limit) {
        // close to right boundary
        // use linear approximation
        lode = 2. + high_slope * (lode_a - 1.);
        dlode = high_slope;
    }
    else {
        const auto lode_b = std::acos(lode_a) / 3.; // theta
        lode = 2. * std::cos(lode_b);               // 2*cos(theta)
        dlode = 2. / 3. * std::sin(lode_b) / std::sqrt((1. - lode_a) * (1. + lode_a));
    }

    const auto square_lode = lode * lode;
    const rowvec dlde = dlode * double_shear * (tensor::stress::lode_der(dev_stress) % tensor::stress::norm_weight).t() * unit_dev_tensor;

    auto ini_f = 0.;
    auto gamma = 0., s = trial_s, p = trial_p;

    mat44 jacobian(fill::none);
    jacobian(0, 0) = 0.;

    vec4 residual, incre;

    mat::fixed<4, 6> left(fill::zeros);

    vec18 data(fill::zeros);
    const auto& f = data(0);
    const auto& pfps = data(1);
    const auto& pfpp = data(2);
    const auto& pfpkp = data(3);
    const auto& gs = data(4);
    const auto& gp = data(5);
    const auto& gg = data(6);
    const auto& pgsps = data(7);
    const auto& pgspp = data(8);
    const auto& pgspkp = data(9);
    const auto& pgpps = data(10);
    const auto& pgppp = data(11);
    const auto& pgppkp = data(12);
    const auto& xh = data(13);
    const auto& dxhdp = data(14);
    const auto& pfpl = data(15);
    // const auto& r = data(16);
    // const auto& drdl = data(17);

    auto counter = 0u;
    auto ref_error = 1.;
    auto try_bisection = false;
    while(true) {
        if(max_iteration == ++counter) {
            if(try_bisection) {
                suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
                return SUANPAN_FAIL;
            }

            try_bisection = true;
            counter = 2u; // avoid initial elastic check

            const auto approx_update = [&](const double gm) {
                gamma = gm;
                s = trial_s - double_shear * gamma * gs;
                p = trial_p - bulk * gamma * gp;
                kp = current_kp + gamma * gg * square_lode / xh;
                compute_plasticity(lode, s, p, kp, data);
                return f;
            };

            auto x1 = 0., f1 = ini_f;
            gamma = f1 / elastic_modulus / elastic_modulus;
            // find a proper bracket
            while(approx_update(gamma) >= 0.) {
                x1 = gamma;
                f1 = f;
                gamma *= 2.;
            }

            ridders(approx_update, x1, f1, gamma, f, tolerance);
        }

        compute_plasticity(lode, s, p, kp, data);

        if(!data.is_finite()) {
            suanpan_error("Non-finite value detected.\n");
            return SUANPAN_FAIL;
        }

        if(1u == counter) {
            if(f < 0.) break;
            ini_f = f;
        }

        residual(0) = f;
        residual(1) = s + double_shear * gamma * gs - trial_s;
        residual(2) = p + bulk * gamma * gp - trial_p;
        residual(3) = xh * (current_kp - kp) + gamma * gg * square_lode;

        jacobian(0, 1) = pfps;
        jacobian(0, 2) = pfpp;
        jacobian(0, 3) = pfpkp;

        jacobian(1, 0) = double_shear * gs;
        jacobian(1, 1) = double_shear * gamma * pgsps + 1.;
        jacobian(1, 2) = double_shear * gamma * pgspp;
        jacobian(1, 3) = double_shear * gamma * pgspkp;

        jacobian(2, 0) = bulk * gp;
        jacobian(2, 1) = bulk * gamma * pgpps;
        jacobian(2, 2) = bulk * gamma * pgppp + 1.;
        jacobian(2, 3) = bulk * gamma * pgppkp;

        jacobian(3, 0) = gg * square_lode;
        jacobian(3, 1) = gamma * square_lode / gg * (gs * pgsps + gp / 3. * pgpps);
        jacobian(3, 2) = gamma * square_lode / gg * (gs * pgspp + gp / 3. * pgppp) + (current_kp - kp) * dxhdp;
        jacobian(3, 3) = gamma * square_lode / gg * (gs * pgspkp + gp / 3. * pgppkp) - xh;

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate + solve_opts::refine)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local plasticity iteration error: {:.5E}.\n", error);

        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) {
            const vec unit_n = n % tensor::stress::norm_weight;

            plastic_strain += gamma * gs * unit_n + gamma * gp / 3. * tensor::unit_tensor2;

            trial_stress = s * n + p * tensor::unit_tensor2;

            mat::fixed<4, 6> right(fill::none);

            right.row(0) = -pfpl * dlde;
            right.row(1) = double_shear * unit_n.t() * unit_dev_tensor;
            right.row(2) = bulk * tensor::unit_tensor2.t();
            right.row(3) = -2. * gamma * gg * lode * dlde;

            if(!solve(left, jacobian, right)) return SUANPAN_FAIL;

            trial_stiffness = n * (left.row(1) - s / trial_s * right.row(1)) + s / trial_s * double_shear * unit_dev_tensor;
            trial_stiffness.row(0) += left.row(2);
            trial_stiffness.row(1) += left.row(2);
            trial_stiffness.row(2) += left.row(2);

            break;
        }

        gamma -= incre(0);
        s -= incre(1);
        p -= incre(2);
        kp -= incre(3);
    }

    //
    // damage part
    //

    vec principal_stress;    // 3
    mat principal_direction; // 3x3
    if(!eig_sym(principal_stress, principal_direction, tensor::stress::to_tensor(trial_stress), "std")) return SUANPAN_FAIL;

    std::vector<uword> tp, cp;
    tp.reserve(3);
    cp.reserve(3);
    for(auto I = 0llu; I < 3llu; ++I)
        if(principal_stress(I) > 0.) tp.emplace_back(I);
        else cp.emplace_back(I);

    const uvec t_pattern(tp), c_pattern(cp);

    const auto aca = accu(square(principal_stress(c_pattern)));
    const auto acb = accu(square(principal_stress));
    const auto ac = aca / acb;

    if(SUANPAN_SUCCESS != compute_damage(gamma, s, p, kp, ac, data)) return SUANPAN_FAIL;

    if(const auto &kdt = trial_history(10), &kdc = trial_history(11); DamageType::NODAMAGE == damage_type || (kdt < e0 && kdc < e0)) return SUANPAN_SUCCESS;

    const auto& omegat = trial_history(16);
    const auto& omegac = trial_history(17);
    const rowvec pot(&data(0), 4);
    const rowvec poc(&data(4), 4);
    const auto& potpac = data(8);
    const auto& potpl = data(9);
    const auto& pocpac = data(10);
    const auto& pocpl = data(11);

    rowvec daca = 2. * principal_stress.t();
    daca(t_pattern).fill(0.);
    const rowvec dac = (daca - 2. * ac * principal_stress.t()) / acb;
    const rowvec dacde = dac * transform::compute_jacobian_nominal_to_principal(principal_direction) * trial_stiffness;

    const rowvec potpe = pot * left + potpac * dacde + potpl * dlde;
    const rowvec pocpe = poc * left + pocpac * dacde + pocpl * dlde;

    const auto damage_t = 1. - omegat;
    const auto damage_c = 1. - omegac;

    if(DamageType::ISOTROPIC == damage_type) {
        trial_stiffness *= damage_t * damage_c;
        trial_stiffness -= trial_stress * (damage_t * pocpe + damage_c * potpe);

        trial_stress *= damage_t * damage_c;
    }
    else if(DamageType::ANISOTROPIC == damage_type) {
        const auto get_fraction = [](const vec& stress) {
            const auto compute_fraction = [&stress](const unsigned i, const unsigned j) {
                const auto &a = stress(i), &b = stress(j);

                return suanpan::approx_equal(a, b, 4) ? a + b <= 0. ? 0. : 2. : 2. * (suanpan::ramp(a) - suanpan::ramp(b)) / (a - b);
            };

            return vec{compute_fraction(0, 1), compute_fraction(1, 2), compute_fraction(2, 0)};
        };

        const mat pnn = transform::eigen_to_tensor_base(principal_direction);

        mat tension_projector = pnn.cols(t_pattern) * pnn.cols(t_pattern).t();
        mat tension_derivative = tension_projector + pnn.tail_cols(3) * diagmat(get_fraction(principal_stress)) * pnn.tail_cols(3).t();

        tension_projector.tail_cols(3) *= 2.;
        tension_derivative.tail_cols(3) *= 2.;

        const vec tension_stress = tension_projector * trial_stress;

        trial_stiffness = damage_c * trial_stiffness - trial_stress * pocpe + (omegac - omegat) * tension_derivative * trial_stiffness + tension_stress * (pocpe - potpe);

        trial_stress *= damage_c;
        trial_stress += (omegac - omegat) * tension_stress;
    }

    return SUANPAN_SUCCESS;
}

int CDPM2::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int CDPM2::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int CDPM2::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

std::vector<vec> CDPM2::record(const OutputType P) {
    if(P == OutputType::DT) return {vec{current_history(16)}};
    if(P == OutputType::DC) return {vec{current_history(17)}};

    return Material3D::record(P);
}

void CDPM2::print() {
    suanpan_info("A concrete damage plasticity model based on the CDPM2 model. doi:10.1016/j.ijsolstr.2013.07.008\n");
}
