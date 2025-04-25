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

#include "MultilinearMises1D.h"

double MultilinearMises1D::compute_k(const double p_strain) const {
    for(uword I = 1; I < backbone.n_rows; ++I)
        if(p_strain <= backbone(I, 0)) return backbone(I - 1, 1) + backbone(I - 1, 2) * (p_strain - backbone(I - 1, 0));
    return backbone(backbone.n_rows - 1, 1);
}

double MultilinearMises1D::compute_dk(const double p_strain) const {
    for(uword I = 1; I < backbone.n_rows; ++I)
        if(p_strain <= backbone(I, 0)) return backbone(I - 1, 2);
    return 0.;
}

double MultilinearMises1D::compute_h(const double) const { return 0.; }

double MultilinearMises1D::compute_dh(const double) const { return 0.; }

MultilinearMises1D::MultilinearMises1D(const unsigned T, const double E, mat&& H, const double R)
    : DataMultilinearMises1D{}
    , NonlinearMises1D(T, E, R) {
    if(H(0, 0) != 0. || H.n_cols != 2) throw std::invalid_argument("first strain should be zero and there should be exact two columns");

    H.resize(H.n_rows, 3);
    H(H.n_rows - 1, 2) = 0.;

    for(uword I = 0; I < H.n_rows - 1; ++I) H(I, 2) = (H(I + 1, 1) - H(I, 1)) / (H(I + 1, 0) - H(I, 0));

    access::rw(backbone) = std::move(H);
}

unique_ptr<Material> MultilinearMises1D::get_copy() { return make_unique<MultilinearMises1D>(*this); }

void MultilinearMises1D::print() {
    suanpan_info("A uniaxial multilinear hardening material using J2 plasticity and associated flow rule.\n");
    NonlinearMises1D::print();
}
