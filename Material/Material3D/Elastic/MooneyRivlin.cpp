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

#include "MooneyRivlin.h"
#include <Toolbox/tensorToolbox.h>

const vec MooneyRivlin::weight{2., 2., 2., 1., 1., 1.};
const vec MooneyRivlin::I1E{2., 2., 2., 0., 0., 0.};
const mat MooneyRivlin::I2EE;
constexpr double MooneyRivlin::one_three = 1. / 3.;
constexpr double MooneyRivlin::two_three = 2. * one_three;
constexpr double MooneyRivlin::four_three = 2. * two_three;
constexpr double MooneyRivlin::five_three = 5. * one_three;
constexpr double MooneyRivlin::eight_nine = two_three * four_three;

MooneyRivlin::MooneyRivlin(const unsigned T, const double KK, const double AA, const double AB, const double R)
    : DataMooneyRivlin{fabs(KK), fabs(AA), fabs(AB)}
    , Material3D(T, R) {}

int MooneyRivlin::initialize(const shared_ptr<DomainBase>&) {
    if(I2EE.is_empty()) {
        mat TI2EE(6, 6, fill::zeros);
        TI2EE(span(0, 2), span(0, 2)).fill(4.);
        TI2EE(0, 0) = TI2EE(1, 1) = TI2EE(2, 2) = 0.;
        TI2EE(3, 3) = TI2EE(4, 4) = TI2EE(5, 5) = -2.;
        access::rw(I2EE) = TI2EE;
    }

    update_trial_status(zeros(6));
    current_stiffness = initial_stiffness = trial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> MooneyRivlin::get_copy() { return make_unique<MooneyRivlin>(*this); }

// takes green strain as input
int MooneyRivlin::update_trial_status(const vec& t_strain) {
    const vec G = weight % (trial_strain = t_strain) + tensor::unit_tensor2;

    const auto &C1 = G(0), &C2 = G(1), &C3 = G(2), &C4 = G(3), &C5 = G(4), &C6 = G(5);

    const auto I1 = C1 + C2 + C3;
    const auto I2 = C1 * C2 + C1 * C3 + C2 * C3 - C4 * C4 - C5 * C5 - C6 * C6;
    const auto I3 = std::max(datum::eps, C1 * C2 * C3 + 2. * C4 * C5 * C6 - C1 * C5 * C5 - C2 * C6 * C6 - C3 * C4 * C4);

    const auto J3M1 = sqrt(I3) - 1.;

    const vec I2E{C2 + C3, C3 + C1, C1 + C2, -C4, -C5, -C6};
    const vec I3E{C2 * C3 - C5 * C5, C3 * C1 - C6 * C6, C1 * C2 - C4 * C4, C5 * C6 - C3 * C4, C6 * C4 - C1 * C5, C4 * C5 - C2 * C6};

    auto W1 = pow(I3, -one_three);
    auto W2 = two_three * I1 * pow(I3, -four_three);
    auto W3 = 2. * W1 * W1;
    auto W4 = four_three * I2 * pow(I3, -five_three);
    auto W5 = pow(I3, -.5);

    const vec J1E = W1 * I1E - W2 * I3E;
    const vec J2E = W3 * I2E - W4 * I3E;
    const vec J3E = W5 * I3E;

    trial_stress = A10 * J1E + A01 * J2E + K * J3M1 * J3E;

    mat I3EE(6, 6, fill::zeros);
    I3EE(1, 2) = I3EE(2, 1) = -2. * (I3EE(4, 4) = -2. * C1);
    I3EE(0, 2) = I3EE(2, 0) = -2. * (I3EE(5, 5) = -2. * C2);
    I3EE(0, 1) = I3EE(1, 0) = -2. * (I3EE(3, 3) = -2. * C3);
    I3EE(2, 3) = I3EE(3, 2) = -2. * (I3EE(4, 5) = I3EE(5, 4) = 2. * C4);
    I3EE(0, 4) = I3EE(4, 0) = -2. * (I3EE(3, 5) = I3EE(5, 3) = 2. * C5);
    I3EE(1, 5) = I3EE(5, 1) = -2. * (I3EE(3, 4) = I3EE(4, 3) = 2. * C6);

    const auto W8 = W5;
    const auto W6 = .5 * W3;
    W1 = two_three * W8;
    W2 *= four_three;
    W3 = .375 * W2;
    W5 = two_three * W4;
    W4 = four_three * W8;
    const auto W7 = .75 * W5;
    const auto W9 = .5 * W8;

    const mat TA = A10 * W1 * J1E + A01 * W4 * J2E;
    const mat TB = TA * J3E.t();

    trial_stiffness = (A10 * W2 + A01 * W5 + K - K * J3M1 * W8) * J3E * J3E.t() + (K * J3M1 * W9 - A10 * W3 - A01 * W7) * I3EE + A01 * W6 * I2EE - TB - TB.t();

    return SUANPAN_SUCCESS;
}

int MooneyRivlin::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    trial_strain.zeros();
    trial_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return SUANPAN_SUCCESS;
}

int MooneyRivlin::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int MooneyRivlin::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void MooneyRivlin::print() { suanpan_info("Mooney-Rivlin material.\n"); }
