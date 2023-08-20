/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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
#include <Toolbox/tensor.h>
#include <Toolbox/utility.h>

const double CDPM2::sqrt_six = std::sqrt(6.);
const double CDPM2::sqrt_three_two = std::sqrt(1.5);
const mat CDPM2::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

void CDPM2::compute_plasticity(const double s, const double p, const double kp, podarray<double>& data) const {
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

    auto qh1 = 1., qh2 = 1.;
    auto dqh1dkp = 0., dqh2dkp = 0.;

    if(kp < 1.) {
        qh1 = qh0 + (1. - qh0) * kp * (kp * (kp - 3.) + 3.) - hp * kp * (kp - 1.) * (kp - 2.);
        dqh1dkp = (3. - 3. * qh0) * pow(kp - 1., 2.) - hp * (3. * kp * (kp - 2.) + 2.);
    }
    else {
        qh2 = 1. + hp * (kp - 1.);
        dqh2dkp = hp;
    }

    const auto ag = 3. * ftfc * qh2 + .5 * m0;
    const auto dagdkp = 3. * ftfc * dqh2dkp;
    const auto cg = qh2 * (1. + ftfc) / 3.;
    const auto dg = log(ag) - log(3. * qh2 + .5 * m0) + lndf;
    const auto dcgdkp = dqh2dkp * (1. + ftfc) / 3.;
    const auto ddgdkp = dagdkp / ag - dqh2dkp / (qh2 + m0 / 6.);
    const auto bg = cg / dg;
    const auto dbgdkp = (dcgdkp - bg * ddgdkp) / dg;

    const auto eg = (p / fc - qh2 * ftfc / 3.) / bg;
    const auto pegpkp = (-ftfc / 3. * dqh2dkp - eg * dbgdkp) / bg;
    const auto pmgpp = ag * exp(eg);

    const auto g3 = (s / sqrt_six + p) / fc;
    const auto g1 = (1. - qh1) * g3 * g3 + sqrt_three_two * s / fc;

    const auto pg3pp = 1. / fc;
    const auto pg3ps = pg3pp / sqrt_six;

    const auto pg2pp = pmgpp / fc;
    const auto pg2ps = m0 * pg3ps;

    const auto pg1pp = (2. - 2. * qh1) * g3 * pg3pp;
    const auto pg1ps = (2. - 2. * qh1) * g3 * pg3ps + sqrt_three_two / fc;
    const auto pg1pkp = -dqh1dkp * g3 * g3;

    f = g1 * g1 + m0 * qh1 * qh1 * qh2 * g3 - qh1 * qh1 * qh2 * qh2;

    pfpp = 2. * g1 * pg1pp + m0 * qh1 * qh1 * qh2 * pg3pp;
    pfps = 2. * g1 * pg1ps + m0 * qh1 * qh1 * qh2 * pg3ps;
    pfpkp = 2. * g1 * pg1pkp + 2. * qh1 * qh2 * (m0 * g3 * dqh1dkp - qh1 * dqh2dkp - qh2 * dqh1dkp) + m0 * qh1 * qh1 * g3 * dqh2dkp;

    gp = 2. * g1 * pg1pp + qh1 * qh1 * pg2pp;
    gs = 2. * g1 * pg1ps + qh1 * qh1 * pg2ps;

    gg = sqrt(gs * gs + gp * gp / 3.);

    pgppp = (4. - 4. * qh1) * (pg1pp * g3 + g1 * pg3pp) + qh1 * qh1 * ag * exp(eg) / bg / fc;
    pgpps = (4. - 4. * qh1) * (pg1ps * g3 + g1 * pg3ps);
    pgppkp = 4. * g3 * (pg1pkp - qh1 * pg1pkp - dqh1dkp * g1) + (2. * dqh1dkp * ag + qh1 * (dagdkp + ag * pegpkp)) * qh1 * exp(eg);

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
        xh = (bh - ah) * exp(-rh / ch);
        dxhdp = xh / ch / fc;
        xh += ah;
    }
    else {
        xh = eh * exp(rh / fh);
        dxhdp = -xh / fh / fc;
        xh += dh;
    }
}

int CDPM2::compute_damage(const double gamma, const double s, const double p, const double kp, const double ac, podarray<double>& data) {
    const auto& gs = data(4);
    const auto& gp = data(5);
    const auto& gg = data(6);
    const auto& pgsps = data(7);
    const auto& pgspp = data(8);
    const auto& pgspkp = data(9);
    const auto& pgpps = data(10);
    const auto& pgppp = data(11);
    const auto& pgppkp = data(12);

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
    const auto ptaps = ptapp / sqrt_six;
    const auto ptbps = sqrt_three_two * e0 / fc;
    const auto term_a = ptaps * s + ptapp * p;
    const auto term_b = ptbps * s;
    const auto term_c = sqrt(term_a * term_a + term_b * term_b);
    ee = term_a + term_c;
    const auto incre_ee = ee - current_ee;
    const auto peeps = (ptaps * ee + term_b * ptbps) / term_c;
    const auto peepp = ptapp * ee / term_c;

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

    // kdt
    auto incre_kdt = 0., pkdtps = 0., pkdtpp = 0.;
    if((et = current_et + incre_ee) > kdt) {
        incre_kdt = et - kdt;
        kdt = et;
        pkdtps = peeps;
        pkdtpp = peepp;
    }

    // kdt1
    auto pkdt1pg = 0., pkdt1ps = 0., pkdt1pp = 0., pkdt1pkp = 0.;
    if(incre_kdt > 0. && kdt > e0) {
        const auto incre_kdt1 = ep / xs;
        kdt1 = current_kdt1 + incre_kdt1;
        pkdt1pg = peppg / xs;
        pkdt1ps = (pepps - incre_kdt1 * pxsps) / xs;
        pkdt1pp = (peppp - incre_kdt1 * pxspp) / xs;
        pkdt1pkp = peppkp / xs;
    }

    // kdt2
    const auto incre_kdt2 = incre_kdt / xs;
    kdt2 = current_kdt2 + incre_kdt2;
    const auto pkdt2ps = (pkdtps - incre_kdt2 * pxsps) / xs;
    const auto pkdt2pp = (pkdtpp - incre_kdt2 * pxspp) / xs;

    // kdc
    auto incre_kdc = 0., pkdcps = 0., pkdcpp = 0., pkdcpac = 0.;
    if((ec = current_ec + ac * incre_ee) > kdc) {
        incre_kdc = ec - kdc;
        kdc = ec;
        pkdcps = ac * peeps;
        pkdcpp = ac * peepp;
        pkdcpac = incre_ee;
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

    podarray<double> datad(3);

    if(SUANPAN_SUCCESS != compute_damage_factor(kdt, kdt1, kdt2, eft, omegat, datad)) return SUANPAN_FAIL;
    const auto& potpkdt = datad(0);
    const auto& potpkdt1 = datad(1);
    const auto& potpkdt2 = datad(2);

    auto& potpg = data(0);
    auto& potps = data(1);
    auto& potpq = data(2);
    auto& potpkp = data(3);

    potpg = potpkdt1 * pkdt1pg;
    potps = potpkdt * pkdtps + potpkdt1 * pkdt1ps + potpkdt2 * pkdt2ps;
    potpq = potpkdt * pkdtpp + potpkdt1 * pkdt1pp + potpkdt2 * pkdt2pp;
    potpkp = potpkdt1 * pkdt1pkp;

    if(SUANPAN_SUCCESS != compute_damage_factor(kdc, kdc1, kdc2, efc, omegac, datad)) return SUANPAN_FAIL;
    const auto& pocpkdc = datad(0);
    const auto& pocpkdc1 = datad(1);
    const auto& pocpkdc2 = datad(2);

    auto& pocpg = data(4);
    auto& pocps = data(5);
    auto& pocpq = data(6);
    auto& pocpkp = data(7);
    auto& pocpac = data(8);

    pocpg = pocpkdc1 * pkdc1pg;
    pocps = pocpkdc * pkdcps + pocpkdc1 * pkdc1ps + pocpkdc2 * pkdc2ps;
    pocpq = pocpkdc * pkdcpp + pocpkdc1 * pkdc1pp + pocpkdc2 * pkdc2pp;
    pocpkp = pocpkdc1 * pkdc1pkp;
    pocpac = pocpkdc * pkdcpac + pocpkdc1 * pkdc1pac + pocpkdc2 * pkdc2pac;

    return SUANPAN_SUCCESS;
}

int CDPM2::compute_damage_factor(const double kd, const double kd1, const double kd2, const double ef, double& omega, podarray<double>& data) const {
    auto& popkd = data(0);
    auto& popkd1 = data(1);
    auto& popkd2 = data(2);

    if(kd < e0) {
        popkd = 0.;
        popkd1 = 0.;
        popkd2 = 0.;
        return SUANPAN_SUCCESS;
    }

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) return SUANPAN_FAIL;

        const auto term_a = ft * exp(-(kd1 + omega * kd2) / ef);
        const auto term_b = (1. - omega) * elastic_modulus;
        const auto jacobian = elastic_modulus * kd - kd2 / ef * term_a;
        const auto incre = (term_a - term_b * kd) / jacobian;

        const auto error = fabs(incre);
        suanpan_debug("Local damage iteration error: {:.5E}.\n", error);

        if(error <= tolerance) {
            popkd = term_b / jacobian;
            popkd1 = term_a / ef / jacobian;
            popkd2 = popkd1 * omega;
            return SUANPAN_SUCCESS;
        }

        omega -= incre;
    }
}

CDPM2::CDPM2(const unsigned T, const double E, const double V, const double FT, const double FC, const double QH0, const double HP, const double DF, const double AH, const double BH, const double CH, const double DH, const double AS, const double EFT, const double EFC, const DamageType DT, const double R)
    : DataCDPM2{fabs(E), fabs(V), fabs(FT), fabs(FC), QH0, HP, DF, AH, BH, CH, DH, AS, fabs(EFT), fabs(EFC)}
    , Material3D(T, R)
    , damage_type(DT) { access::rw(tolerance) = 1E-13; }

int CDPM2::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(18);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> CDPM2::get_copy() { return make_unique<CDPM2>(*this); }

double CDPM2::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return elastic_modulus;
    if(ParameterType::SHEARMODULUS == P || ParameterType::G == P) return .5 * double_shear;
    if(ParameterType::BULKMODULUS == P) return elastic_modulus / (3. - 6. * poissons_ratio);
    if(ParameterType::POISSONSRATIO == P) return poissons_ratio;
    return 0.;
}

int CDPM2::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= tolerance) return SUANPAN_SUCCESS;

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

    auto gamma = 0., s = trial_s, p = trial_p;

    mat jacobian(4, 4, fill::none), left(4, 6, fill::zeros);
    jacobian(0, 0) = 0.;

    vec residual(4), incre;

    podarray<double> data(15);
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

    auto counter = 0u;

    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        compute_plasticity(s, p, kp, data);

        if(1 == counter && f < 0.) break;

        residual(0) = f;
        residual(1) = s + double_shear * gamma * gs - trial_s;
        residual(2) = p + bulk * gamma * gp - trial_p;
        residual(3) = xh * (current_kp - kp) + gamma * gg;

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

        jacobian(3, 0) = gg;
        jacobian(3, 1) = gamma / gg * (gs * pgsps + gp / 3. * pgpps);
        jacobian(3, 2) = gamma / gg * (gs * pgspp + gp / 3. * pgppp) + (current_kp - kp) * dxhdp;
        jacobian(3, 3) = gamma / gg * (gs * pgspkp + gp / 3. * pgppkp) - xh;

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = norm(residual);
        suanpan_debug("Local plasticity iteration error: {:.5E}.\n", error);

        if(error <= tolerance) {
            const vec unit_n = n % tensor::stress::norm_weight;

            plastic_strain += gamma * gs * unit_n + gamma * gp / 3. * tensor::unit_tensor2;

            trial_stress = s * n + p * tensor::unit_tensor2;

            mat right(4, 6, fill::none);

            right.row(0).zeros();
            right.row(1) = double_shear * unit_n.t() * unit_dev_tensor;
            right.row(2) = bulk * tensor::unit_tensor2.t();
            right.row(3).zeros();

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

    vector<uword> tp, cp;
    tp.reserve(3);
    cp.reserve(3);
    for(auto I = 0llu; I < 3llu; ++I)
        if(principal_stress(I) > 0.) tp.emplace_back(I);
        else cp.emplace_back(I);

    const uvec t_pattern(tp), c_pattern(cp);

    const auto aca = accu(square(principal_stress(c_pattern)));
    const auto acb = accu(square(principal_stress));
    const auto ac = aca / acb;
    rowvec daca = 2. * principal_stress.t();
    daca(t_pattern).fill(0.);
    const rowvec dac = (daca - 2. * ac * principal_stress.t()) / acb;

    if(SUANPAN_SUCCESS != compute_damage(gamma, s, p, kp, ac, data)) return SUANPAN_FAIL;

    if(DamageType::NODAMAGE == damage_type) return SUANPAN_SUCCESS;

    const auto& omegat = trial_history(16);
    const auto& omegac = trial_history(17);
    const rowvec pot(&data(0), 4);
    const rowvec poc(&data(4), 4);
    const auto& pocpac = data(8);

    const rowvec potpe = pot * left;
    const rowvec pocpe = poc * left + pocpac * dac * transform::compute_jacobian_nominal_to_principal(principal_direction) * trial_stiffness;

    if(DamageType::ISOTROPIC == damage_type) {
        trial_stiffness *= (1. - omegat) * (1. - omegac);
        trial_stiffness -= trial_stress * ((1. - omegat) * pocpe + (1. - omegac) * potpe);

        trial_stress *= (1. - omegat) * (1. - omegac);
    }
    else if(DamageType::ANISOTROPIC == damage_type) {
        const auto get_fraction = [](const vec& p_stress) {
            const auto compute_fraction = [](const double a, const double b) { return suanpan::approx_equal(a, b, 4) ? a + b <= 0. ? 0. : 2. : 2. * (suanpan::ramp(a) - suanpan::ramp(b)) / (a - b); };

            return vec{compute_fraction(p_stress(0), p_stress(1)), compute_fraction(p_stress(1), p_stress(2)), compute_fraction(p_stress(2), p_stress(0))};
        };

        const mat pnn = transform::eigen_to_tensor_base(principal_direction);

        mat tension_projector = pnn.cols(t_pattern) * pnn.cols(t_pattern).t();
        mat tension_derivative = tension_projector + pnn.tail_cols(3) * diagmat(get_fraction(principal_stress)) * pnn.tail_cols(3).t();

        tension_projector.tail_cols(3) *= 2.;
        tension_derivative.tail_cols(3) *= 2.;

        const vec tension_stress = tension_projector * trial_stress;

        trial_stiffness = (1. - omegac) * trial_stiffness - trial_stress * pocpe + (omegac - omegat) * tension_derivative * trial_stiffness + tension_stress * (pocpe - potpe);

        trial_stress *= 1. - omegac;
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

vector<vec> CDPM2::record(const OutputType T) {
    if(T == OutputType::KAPPAP) return {vec{current_history(0)}};
    if(T == OutputType::DT) return {vec{current_history(16)}};
    if(T == OutputType::DC) return {vec{current_history(17)}};

    return Material3D::record(T);
}

void CDPM2::print() {
    suanpan_info("A concrete damage plasticity model based on the CDPM2 model. doi: 10.1016/j.ijsolstr.2013.07.008\n");
}
