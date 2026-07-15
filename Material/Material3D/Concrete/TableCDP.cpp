/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#include "TableCDP.h"

#include <ranges>

pod2 TableCDP::interpolate(const double kappa, const mat& table) {
    const auto indices = std::views::iota(uword{0}, table.n_rows);

    const auto right = *std::ranges::upper_bound(indices, kappa, {}, [&](const uword idx) { return table(idx, 0); }), left = right - 1;

    const double x0 = table(left, 0), x1 = table(right, 0), y0 = table(left, 1), y1 = table(right, 1);

    const auto rate = (y1 - y0) / (x1 - x0);

    return {rate * (kappa - x0) + y0, rate};
}

pod6 TableCDP::compute_backbone(const double kappa, const double sign, const mat& stress_table, const mat& damage_table) {
    pod6 out;

    if(kappa < 0.) {
        out[3] = 0.; // \md{d}
        out[0] = 0.; // d
    }
    else if(kappa >= 1.) {
        out[3] = 0.;              // \md{d}
        out[0] = 1. - datum::eps; // d
    }
    else if(const auto last_row = damage_table.n_rows - 1; kappa >= damage_table(last_row, 0)) {
        // connecting the last point in the table to (1,1)
        out[3] = (1. - damage_table(last_row, 1)) / (1. - damage_table(last_row, 0)); // \md{d}
        out[0] = std::min(1. - datum::eps, 1. + (kappa - 1.) * out[3]);               // d
    }
    else {
        const auto damage = interpolate(kappa, damage_table);
        out[0] = damage[0]; // d
        out[3] = damage[1]; // \md{d}
    }

    if(kappa < 0.) {
        out[1] = stress_table(0, 0); // f
        out[4] = 0.;                 // \md{f}
    }
    else if(kappa >= 1.) {
        out[1] = 0.; // f
        out[4] = 0.; // \md{f}
    }
    else if(const auto last_row = stress_table.n_rows - 1; kappa >= stress_table(last_row, 0)) {
        // connecting the last point in the table to (1,0)
        out[4] = stress_table(last_row, 1) / (stress_table(last_row, 0) - 1.); // \md{f}
        out[1] = out[4] * (kappa - 1.);                                        // f
    }
    else {
        const auto stress = interpolate(kappa, stress_table);
        out[1] = stress[0]; // f
        out[4] = stress[1]; // \md{f}
    }

    out[1] *= sign;
    out[4] *= sign;

    out[2] = out[1] / (1. - out[0]);                                                  // \bar{f}
    out[5] = ((1. - out[0]) * out[4] + out[1] * out[3]) * std::pow(1. - out[0], -2.); // \md{\bar{f}}

    return out;
}

pod6 TableCDP::compute_tension_backbone(const double kappa) const { return compute_backbone(kappa, 1., t_table, dt_table); }

pod6 TableCDP::compute_compression_backbone(const double kappa) const { return compute_backbone(kappa, -1., c_table, dc_table); }

TableCDP::TableCDP(const unsigned T, const double E, const double V, mat&& TT, mat&& CT, mat&& TDT, mat&& CDT, const double GT, const double GC, const double AP, const double BC, const double S, const double R)
    : NonlinearCDP(T, E, V, GT, GC, AP, BC, S, R)
    , t_table(std::move(TT))
    , c_table(std::move(CT))
    , dt_table(std::move(TDT))
    , dc_table(std::move(CDT)) {}

unique_ptr<Material> TableCDP::unique_copy() { return std::make_unique<TableCDP>(*this); }
