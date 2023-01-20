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

#include "MultilinearElastic1D.h"

MultilinearElastic1D::MultilinearElastic1D(const unsigned T, mat&& H, const double R)
    : Material1D(T, R) {
    if(H.n_cols != 2) throw invalid_argument("there should be exact two columns.\n");
    H.resize(H.n_rows, 3);
    H(0, 2) = H(0, 1) / H(0, 0);
    for(unsigned I = 1; I < H.n_rows; ++I) H(I, 2) = (H(I, 1) - H(I - 1, 1)) / (H(I, 0) - H(I - 1, 0));
    access::rw(backbone) = std::forward<mat>(H);
}

int MultilinearElastic1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = backbone(0, 1) / backbone(0, 0);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> MultilinearElastic1D::get_copy() { return make_unique<MultilinearElastic1D>(*this); }

int MultilinearElastic1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    const auto abs_strain = fabs(trial_strain(0));

    auto flag = true;
    for(unsigned I = 0; I < backbone.n_rows; ++I)
        if(abs_strain < backbone(I, 0)) {
            trial_stress = backbone(I, 1) + (abs_strain - backbone(I, 0)) * (trial_stiffness = backbone(I, 2));
            flag = false;
            break;
        }

    if(flag) {
        trial_stiffness = 0.;
        trial_stress = backbone(backbone.n_rows - 1, 1);
    }

    if(trial_strain(0) < 0.) trial_stress = -trial_stress;

    return SUANPAN_SUCCESS;
}

int MultilinearElastic1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_stiffness = initial_stiffness;
    return reset_status();
}

int MultilinearElastic1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int MultilinearElastic1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void MultilinearElastic1D::print() {
    suanpan_info("A multilinear elastic material model.\n");
    Material1D::print();
}
