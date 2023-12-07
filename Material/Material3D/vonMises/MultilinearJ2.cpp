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

#include "MultilinearJ2.h"

double MultilinearJ2::compute_k(const double p_strain) const {
    for(unsigned I = 1; I < backbone.n_rows; ++I) if(p_strain <= backbone(I, 0)) return backbone(I - 1, 1) + backbone(I - 1, 2) * (p_strain - backbone(I - 1, 0));
    return backbone(backbone.n_rows - 1, 1);
}

double MultilinearJ2::compute_dk(const double p_strain) const {
    for(unsigned I = 1; I < backbone.n_rows; ++I) if(p_strain <= backbone(I, 0)) return backbone(I - 1, 2);
    return 0.;
}

double MultilinearJ2::compute_h(const double) const { return 0.; }

double MultilinearJ2::compute_dh(const double) const { return 0.; }

MultilinearJ2::MultilinearJ2(const unsigned T, const double E, const double V, mat&& H, const double R)
    : NonlinearJ2(T, E, V, R) {
    if(H(0, 0) != 0. || H.n_cols != 2) throw invalid_argument("first strain should be zero and there should be exact two columns");

    H.resize(H.n_rows, 3);
    H(H.n_rows - 1, 2) = 0.;

    for(unsigned I = 0; I < H.n_rows - 1; ++I) H(I, 2) = (H(I + 1llu, 1) - H(I, 1)) / (H(I + 1llu, 0) - H(I, 0));

    access::rw(backbone) = std::move(H);
}

unique_ptr<Material> MultilinearJ2::get_copy() { return make_unique<MultilinearJ2>(*this); }

void MultilinearJ2::print() {
    suanpan_info("A 3D multilinear hardening model.\nE = {:.4E}\t\\sigma_y = {:.4E}.\n", elastic_modulus, backbone(0, 1));
}
