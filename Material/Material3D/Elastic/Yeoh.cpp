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

#include "Yeoh.h"
#include <Toolbox/tensor.h>

const vec Yeoh::weight{2., 2., 2., 1., 1., 1.};
const vec Yeoh::I1E{2., 2., 2., 0., 0., 0.};

vec Yeoh::compute_derivative(const double J1M3, const double J3M1) const {
    const auto J3M12 = J3M1 * J3M1;

    vec D(4, fill::zeros);

    auto &DWDJ1 = D(0), &DWDJ3 = D(1), &DDWDDJ1 = D(2), &DDWDDJ3 = D(3);

    auto TMP = 1., IDX = 1.;
    for(auto X = 0llu; X < A0.n_elem; ++X) {
        DWDJ1 += A0(X) * IDX * TMP;
        ++IDX;
        TMP *= J1M3;
    }

    TMP = 1.;
    IDX = 1.;
    for(auto X = 1llu; X < A0.n_elem; ++X) {
        DDWDDJ1 += A0(X) * IDX * (IDX + 1.) * TMP;
        ++IDX;
        TMP *= J1M3;
    }

    TMP = J3M1;
    IDX = 2.;
    for(auto X = 0llu; X < A1.n_elem; ++X) {
        DWDJ3 += A1(X) * IDX * TMP;
        TMP *= J3M12;
        IDX += 2.;
    }

    TMP = 1.;
    IDX = 1.;
    for(auto X = 0llu; X < A1.n_elem; ++X) {
        DDWDDJ3 += A1(X) * IDX * (IDX + 1.) * TMP;
        IDX += 2.;
        TMP *= J3M12;
    }

    return D;
}

Yeoh::Yeoh(const unsigned T, vec&& CC, vec&& KK, const double R)
    : DataYeoh{std::move(CC), std::move(KK)}
    , Material3D(T, R) {}

int Yeoh::initialize(const shared_ptr<DomainBase>&) {
    update_trial_status(zeros(6));
    current_stiffness = initial_stiffness = trial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Yeoh::get_copy() { return make_unique<Yeoh>(*this); }

// takes green strain as input
int Yeoh::update_trial_status(const vec& t_strain) {
    const vec G = weight % (trial_strain = t_strain) + tensor::unit_tensor2;

    const auto &C1 = G(0), &C2 = G(1), &C3 = G(2), &C4 = G(3), &C5 = G(4), &C6 = G(5);

    const auto I1 = C1 + C2 + C3;
    const auto I3 = std::max(datum::eps, C1 * C2 * C3 + 2. * C4 * C5 * C6 - C1 * C5 * C5 - C2 * C6 * C6 - C3 * C4 * C4);

    const vec I3E = 2. * vec{C2 * C3 - C5 * C5, C3 * C1 - C6 * C6, C1 * C2 - C4 * C4, C5 * C6 - C3 * C4, C6 * C4 - C1 * C5, C4 * C5 - C2 * C6};

    const auto W1 = pow(I3, -one_three);
    const auto W2 = one_three * I1 * pow(I3, -four_three);
    const auto W5 = .5 * pow(I3, -.5);

    const vec J1E = W1 * I1E - W2 * I3E;
    const vec J3E = W5 * I3E;

    const auto D = compute_derivative(I1 * W1 - 3., sqrt(I3) - 1.);

    const auto &DWDJ1 = D(0), &DWDJ3 = D(1), &DDWDDJ1 = D(2), &DDWDDJ3 = D(3);

    trial_stress = DWDJ1 * J1E + DWDJ3 * J3E;

    mat I3EE(6, 6, fill::zeros);
    I3EE(1, 2) = I3EE(2, 1) = -2. * (I3EE(4, 4) = -2. * C1);
    I3EE(0, 2) = I3EE(2, 0) = -2. * (I3EE(5, 5) = -2. * C2);
    I3EE(0, 1) = I3EE(1, 0) = -2. * (I3EE(3, 3) = -2. * C3);
    I3EE(2, 3) = I3EE(3, 2) = -2. * (I3EE(4, 5) = I3EE(5, 4) = 2. * C4);
    I3EE(0, 4) = I3EE(4, 0) = -2. * (I3EE(3, 5) = I3EE(5, 3) = 2. * C5);
    I3EE(1, 5) = I3EE(5, 1) = -2. * (I3EE(3, 4) = I3EE(4, 3) = 2. * C6);

    const auto P1 = 2. * four_three * W2 * DWDJ1;
    const auto P2 = four_three * W5 * DWDJ1;
    const auto P3 = W2 * DWDJ1;
    const auto P4 = W5 * DWDJ3;

    trial_stiffness = (P1 + DDWDDJ3 - 2. * P4) * J3E * J3E.t() + (P4 - P3) * I3EE + (DDWDDJ1 * J1E - P2 * J3E) * J1E.t() - P2 * J1E * J3E.t();

    return SUANPAN_SUCCESS;
}

int Yeoh::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    trial_strain.zeros();
    trial_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return SUANPAN_SUCCESS;
}

int Yeoh::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int Yeoh::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void Yeoh::print() {
    suanpan_info("A Yeoh material model.\n");
}
