/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

// ReSharper disable IdentifierTypo
#include "NM3D2.h"

double NM3D2::evaluate(const double p, const double ms, const double mw, const mat& weight) const {
    double value = weight(0);

    if(weight(1) > 0.) value *= pow(p, weight(1));
    if(weight(2) > 0.) value *= pow(ms, weight(2));
    if(weight(3) > 0.) value *= pow(mw, weight(3));

    return value;
}

vec NM3D2::differentiate(const mat& weight, uword location, const uword order) {
    ++location;

    vec weight_out = weight.as_col();

    weight_out(location) = weight(location) - static_cast<double>(order);

    if(weight_out(location) < 0.) weight_out.zeros();
    else for(auto I = static_cast<uword>(weight(location)); I > static_cast<uword>(weight_out(location)); --I) weight_out(0) *= static_cast<double>(I);

    return weight_out;
}

double NM3D2::compute_h(const double alpha) const { return std::max(0., c + h * alpha); }

double NM3D2::compute_dh(const double alpha) const { return (c + h * alpha > 0.) * h; }

double NM3D2::compute_f(const vec& s) const {
    const auto p = s(0) / yield_force(0);
    const auto ms = s(1) / yield_force(1);
    const auto mw = s(2) / yield_force(2);

    auto f = 0.;
    for(auto I = 0llu; I < para_set.n_rows; ++I) f += evaluate(p, ms, mw, para_set.row(I));

    return f;
}

vec NM3D2::compute_df(const vec& s) const {
    const auto p = s(0) / yield_force(0);
    const auto ms = s(1) / yield_force(1);
    const auto mw = s(2) / yield_force(2);

    vec df(3, fill::zeros);

    for(auto I = 0llu; I < para_set.n_rows; ++I) for(auto J = 0llu; J < df.n_elem; ++J) df(J) += evaluate(p, ms, mw, differentiate(para_set.row(I), J, 1));

    return df / yield_force;
}

mat NM3D2::compute_ddf(const vec& s) const {
    const auto p = s(0) / yield_force(0);
    const auto ms = s(1) / yield_force(1);
    const auto mw = s(2) / yield_force(2);

    mat ddf(3, 3, fill::zeros);

    for(auto I = 0llu; I < para_set.n_rows; ++I)
        for(auto J = 0llu; J < ddf.n_rows; ++J) {
            const auto dfj = differentiate(para_set.row(I), J, 1);
            for(auto K = 0llu; K < ddf.n_cols; ++K) ddf(J, K) += evaluate(p, ms, mw, differentiate(dfj, K, 1));
        }

    return ddf / (yield_force * yield_force.t());
}

NM3D2::NM3D2(const unsigned T, const double EEA, const double EEIS, const double EEIW, const double NP, const double MSP, const double MWP, const double CC, const double HH, const double KK, const double LD, mat&& PS)
    : NonlinearNM3D(T, EEA, EEIS, EEIW, KK, LD)
    , para_set(PS.empty() ? mat{{1.15, 2., 0., 0.}, {1., 0., 2., 0.}, {1., 0., 0., 4.}, {3.67, 2., 2., 0.}, {3., 6., 0., 2.}, {4.65, 0., 4., 2.}} : std::forward<mat>(PS))
    , yield_force{NP, MSP, MWP}
    , c(CC)
    , h(HH)
    , k(KK) {}

int NM3D2::initialize(const shared_ptr<DomainBase>& D) {
    if(SUANPAN_SUCCESS != NonlinearNM3D::initialize(D)) return SUANPAN_FAIL;

    initialize_weight(yield_force, k);

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> NM3D2::get_copy() { return make_unique<NM3D2>(*this); }

void NM3D2::print() {
    suanpan_info("A N-M based section. doi: 10.1007/978-94-007-6573-3_3\n");
    current_deformation.head(5).t().print("Local Deformation:");
    current_resistance.head(5).t().print("Local Resistance:");
}
