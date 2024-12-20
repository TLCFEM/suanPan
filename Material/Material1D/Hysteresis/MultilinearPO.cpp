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

#include "MultilinearPO.h"

podarray<double> MultilinearPO::compute_tension_initial_reverse() const {
    podarray<double> response(2);

    response(0) = t_backbone(0, 0);
    response(1) = t_backbone(0, 1);

    return response;
}

podarray<double> MultilinearPO::compute_compression_initial_reverse() const {
    podarray<double> response(2);

    response(0) = c_backbone(0, 0);
    response(1) = c_backbone(0, 1);

    return response;
}

podarray<double> MultilinearPO::compute_tension_backbone(const double t_strain) const {
    podarray<double> response(2);

    uword IDX = 0;
    while(IDX < t_backbone.n_rows && t_backbone(IDX, 0) < t_strain) ++IDX;

    if(0 == IDX) {
        response(1) = t_backbone(0, 1) / t_backbone(0, 0);
        response(0) = response(1) * t_strain;
    }
    else if(t_backbone.n_rows == IDX) {
        response(1) = 1E-11;
        response(0) = t_backbone(t_backbone.n_rows - 1, 1);
    }
    else {
        response(1) = (t_backbone(IDX, 1) - t_backbone(IDX - 1, 1)) / (t_backbone(IDX, 0) - t_backbone(IDX - 1, 0));
        response(0) = t_backbone(IDX - 1, 1) + (t_strain - t_backbone(IDX - 1, 0)) * response(1);
    }

    return response;
}

podarray<double> MultilinearPO::compute_compression_backbone(const double t_strain) const {
    podarray<double> response(2);

    uword IDX = 0;
    while(IDX < c_backbone.n_rows && c_backbone(IDX, 0) > t_strain) ++IDX;

    if(0 == IDX) {
        response(1) = c_backbone(0, 1) / c_backbone(0, 0);
        response(0) = response(1) * t_strain;
    }
    else if(c_backbone.n_rows == IDX) {
        response(1) = 0.;
        response(0) = c_backbone(c_backbone.n_rows - 1, 1);
    }
    else {
        response(1) = (c_backbone(IDX, 1) - c_backbone(IDX - 1, 1)) / (c_backbone(IDX, 0) - c_backbone(IDX - 1, 0));
        response(0) = c_backbone(IDX - 1, 1) + (t_strain - c_backbone(IDX - 1, 0)) * response(1);
    }

    return response;
}

MultilinearPO::MultilinearPO(const int T, mat&& TB, mat&& CB, const double R)
    : DataMultilinearPO{std::move(TB), std::move(CB)}
    , PeakOriented(T, R) {}

int MultilinearPO::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = t_backbone(0, 1) / t_backbone(0, 0);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> MultilinearPO::get_copy() { return make_unique<MultilinearPO>(*this); }

void MultilinearPO::print() {
    suanpan_info("A multilinear peak oriented hysteresis model.\n");
    PeakOriented::print();
}
