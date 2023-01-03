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

#include "BWBN.h"

BWBN::BWBN(const unsigned T, vec&& P, const double R)
    : DataBWBN{std::forward<vec>(P)}
    , Material1D(T, R) {}

int BWBN::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(2);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> BWBN::get_copy() { return make_unique<BWBN>(*this); }

int BWBN::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    const auto n_strain = incre_strain(0) / yield_strain;

    trial_history = current_history;
    const auto& current_z = current_history(0); // z
    const auto& current_e = current_history(1); // energy
    auto& z = trial_history(0);                 // z
    auto& e = trial_history(1);                 // energy

    auto incre = .5 * n_strain;
    unsigned counter = 0;
    while(true) {
        if(max_iteration == ++counter) {
            SP_E("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        z += incre;

        e = current_e + .5 * modulus_c * (z + current_z) * n_strain;

        const auto pepz = .5 * modulus_c * n_strain;

        const auto mu = nu_i + nu_rate * e;
        const auto eta = eta_i + eta_rate * e;
        const auto a = 1. - a_rate * e;

        const auto pmupz = nu_rate * pepz;
        const auto petapz = eta_rate * pepz;
        const auto papz = -a_rate * pepz;

        const auto za = (1. - exp(-p * e)) * zeta;
        const auto zb = (phi_i + phi_rate * e) * (lambda + za);
        const auto zm = pow(mu, -1. / n);

        const auto pzapz = zeta * p * exp(-p * e) * pepz;
        const auto pzbpz = (phi_i + phi_rate * e) * pzapz + (lambda + za) * phi_rate * pepz;
        const auto pzmpz = -zm * pmupz / mu / n;

        double f_term, pfpz;
        if(n_strain >= 0.) {
            f_term = (z - q * zm) / zb;
            pfpz = (-f_term * pzbpz + 1. - q * pzmpz) / zb;
        }
        else {
            f_term = (-z - q * zm) / zb;
            pfpz = (-f_term * pzbpz - 1. - q * pzmpz) / zb;
        }

        const auto e_term = exp(-f_term * f_term);

        const auto h = 1. - za * e_term;

        const auto phpz = e_term * (2. * za * f_term * pfpz - pzapz);

        const auto p_term = (z * n_strain >= 0. ? 1. : 1. - 2. * beta) * pow(std::max(datum::eps, fabs(z)), n);
        const auto t_term = mu * n_strain * p_term;

        const auto ptpz = n_strain * (mu * n * p_term / z + p_term * pmupz);

        const auto residual = (z - current_z) * eta + (t_term - a * n_strain) * h;
        const auto jacobian = eta + (z - current_z) * petapz + (ptpz - papz * n_strain) * h + (t_term - a * n_strain) * phpz;

        const auto error = fabs(incre = -residual / jacobian);

        SP_D("Local iteration error: {:.5E}.\n", error);

        if(error <= tolerance) {
            trial_stress = modulus_a * trial_strain + modulus_b * z;

            const auto pepn = .5 * modulus_c * (z + current_z);

            const auto pmupn = nu_rate * pepn;
            const auto petapn = eta_rate * pepn;
            const auto papn = -a_rate * pepn;

            const auto pzapn = zeta * p * exp(-p * e) * pepn;
            const auto pzbpn = (phi_i + phi_rate * e) * pzapn + (lambda + za) * phi_rate * pepn;
            const auto pzmpn = -zm * pmupn / mu / n;

            const auto pfpn = (-f_term * pzbpn - q * pzmpn) / zb;

            const auto phpn = e_term * (2. * za * f_term * pfpn - pzapn);
            const auto ptpn = (pmupn * n_strain + mu) * p_term;

            const auto left = (z - current_z) * petapn + (ptpn - papn * n_strain - a) * h + (t_term - a * n_strain) * phpn;

            trial_stiffness = modulus_a - modulus_b / yield_strain * left / jacobian;

            return SUANPAN_SUCCESS;
        }
    }
}

int BWBN::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int BWBN::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    current_history = trial_history;
    return SUANPAN_SUCCESS;
}

int BWBN::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    trial_history = current_history;
    return SUANPAN_SUCCESS;
}

void BWBN::print() {
    sp_info("A BWBN material model.\n");
    Material1D::print();
}
