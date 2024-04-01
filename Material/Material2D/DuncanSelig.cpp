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
#include <Toolbox/tensor.h>

DuncanSelig::DuncanSelig(const unsigned T, const double E, const double R)
    : Material2D(T, PlaneType::E, R)
    , elastic_modulus(fabs(E)) {}

int DuncanSelig::initialize(const shared_ptr<DomainBase>& D) {
    initial_stiffness.zeros(3, 3);
    initial_stiffness(2, 2) = shear_modulus = .5 * (initial_stiffness(0, 0) = initial_stiffness(1, 1) = elastic_modulus);

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> DuncanSelig::get_copy() { return make_unique<DuncanSelig>(*this); }

int DuncanSelig::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    double ini_phi, ten_fold_phi_diff, p_atm, r_f, cohesion, ref_elastic, n;

    double s1, s3;

    // for elastic modulus

    const auto phi = ini_phi - ten_fold_phi_diff * log10(s3 / p_atm);
    const auto dphids3 = -ten_fold_phi_diff / (s3 * log(10));

    const auto denom = 1. - sin(phi);
    const auto max_dev_stress = 2. / r_f * (cohesion * cos(phi) + s3 * sin(phi)) / denom;
    const auto pmdspphi = 2. / r_f * (s3 * cos(phi) / denom / denom + cohesion / denom);
    const auto pmdsps3 = 2. / r_f * sin(phi) / denom;
    const auto dmdsds3 = pmdspphi * dphids3 + pmdsps3;

    const auto dev_stress = s1 - s3;
    const auto pdsps1 = 1.;
    const auto pdsps3 = -1.;

    const auto ini_elastic = ref_elastic * pow(s3 / p_atm, n);
    const auto deids3 = n * ini_elastic / s3;

    const auto pepei = pow(1. - dev_stress / max_dev_stress, 2.);
    const auto elastic = ini_elastic * pepei;
    const auto pepds = -2. * ini_elastic * (1. - dev_stress / max_dev_stress) / max_dev_stress;
    const auto pepmds = 2. * ini_elastic * (1. - dev_stress / max_dev_stress) * dev_stress / max_dev_stress / max_dev_stress;

    const auto peps1 = pepds * pdsps1;
    const auto peps3 = pepei * deids3 + pepds * pdsps3 + pepmds * dmdsds3;

    // for bulk modulus

    return SUANPAN_SUCCESS;
}

int DuncanSelig::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return SUANPAN_SUCCESS;
}

int DuncanSelig::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int DuncanSelig::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void DuncanSelig::print() {
    suanpan_info("The Duncan-Selig soil model.\n");
    suanpan_info("Strain:", current_strain);
    suanpan_info("Stress:", current_stress);
}
