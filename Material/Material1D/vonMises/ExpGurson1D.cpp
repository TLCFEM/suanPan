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

#include "ExpGurson1D.h"

ExpGurson1D::ExpGurson1D(const unsigned T, const double E, const double V, const double YS, const double N, const double Q1, const double Q2, const double FN, const double SN, const double EN, const double R)
    : DataExpGurson1D{fabs(YS), std::min(1., N)}
    , NonlinearGurson1D(T, E, V, Q1, Q2, FN, SN, EN, R) {}

vec ExpGurson1D::compute_hardening(const double plastic_strain) const {
    auto k = 1.;
    double pow_term;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            k = pow_term = 1.;
            break;
        }

        const auto tmp_term = k + para_c * plastic_strain;
        pow_term = std::pow(tmp_term, n - 1.);
        const auto residual = k - pow_term * tmp_term;
        const auto incre = residual / (1. - n * pow_term);
        const auto error = std::fabs(incre);
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error <= tolerance || (std::fabs(residual) < tolerance && counter > 5u)) break;

        k -= incre;
    }

    pow_term *= n;

    return vec{k, pow_term * para_c / (1. - pow_term)} * yield_stress;
}

unique_ptr<Material> ExpGurson1D::get_copy() { return std::make_unique<ExpGurson1D>(*this); }

void ExpGurson1D::print() {
    suanpan_info("A uniaxial Gurson model.\n");
    NonlinearGurson1D::print();
}
