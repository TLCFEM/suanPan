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

#include "MultilinearOO.h"

pod2 MultilinearOO::compute_tension_backbone(const double t_strain) const {
    pod2 response;

    uword IDX = 0;
    while(IDX < t_backbone.n_rows && t_backbone(IDX, 0) < t_strain) ++IDX;

    if(0 == IDX) {
        response[1] = t_backbone(0, 1) / t_backbone(0, 0);
        response[0] = response[1] * t_strain;
    }
    else if(t_backbone.n_rows == IDX) {
        response[1] = 1E-11;
        response[0] = t_backbone(t_backbone.n_rows - 1, 1);
    }
    else {
        response[1] = (t_backbone(IDX, 1) - t_backbone(IDX - 1, 1)) / (t_backbone(IDX, 0) - t_backbone(IDX - 1, 0));
        response[0] = t_backbone(IDX - 1, 1) + (t_strain - t_backbone(IDX - 1, 0)) * response[1];
    }

    return response;
}

pod2 MultilinearOO::compute_compression_backbone(const double t_strain) const {
    pod2 response;

    uword IDX = 0;
    while(IDX < c_backbone.n_rows && c_backbone(IDX, 0) > t_strain) ++IDX;

    if(0 == IDX) {
        response[1] = c_backbone(0, 1) / c_backbone(0, 0);
        response[0] = response[1] * t_strain;
    }
    else if(c_backbone.n_rows == IDX) {
        response[1] = 0.;
        response[0] = c_backbone(c_backbone.n_rows - 1, 1);
    }
    else {
        response[1] = (c_backbone(IDX, 1) - c_backbone(IDX - 1, 1)) / (c_backbone(IDX, 0) - c_backbone(IDX - 1, 0));
        response[0] = c_backbone(IDX - 1, 1) + (t_strain - c_backbone(IDX - 1, 0)) * response[1];
    }

    return response;
}

MultilinearOO::MultilinearOO(const int T, mat&& TB, mat&& CB, const double R)
    : DataMultilinearOO{std::move(TB), std::move(CB)}
    , OriginOriented(T, R) {}

int MultilinearOO::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = t_backbone(0, 1) / t_backbone(0, 0);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> MultilinearOO::get_copy() { return make_unique<MultilinearOO>(*this); }

void MultilinearOO::print() {
    suanpan_info("A multilinear origin oriented hysteresis model.\n");
    OriginOriented::print();
}
