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

    const auto right = *std::ranges::lower_bound(indices, kappa, {}, [&](const uword idx) { return table(idx, 0); });
    const auto left = right - 1;

    const double x0 = table(left, 0);
    const double x1 = table(right, 0);
    const double y0 = table(left, 1);
    const double y1 = table(right, 1);

    const auto rate = (y1 - y0) / (x1 - x0);
    const auto value = rate * (kappa - x0) + y0;

    return {value, rate};
}

pod6 TableCDP::compute_tension_backbone(const double kappa) const {
    pod6 out;

    if(kappa < dt_table(0, 0)) {
        out[3] = dt_table(0, 1) / dt_table(0, 0); // \md{d}
        out[0] = kappa * out[3];                  // d
    }
    else if(kappa > dt_table(dt_table.n_rows - 1, 0)) {
        out[3] = (1. - dt_table(dt_table.n_rows - 1, 1)) / (1. - dt_table(dt_table.n_rows - 1, 0)); // \md{d}
        out[0] = std::min(1., 1. - (1. - kappa) * out[3]);                                          // d
    }
    else
        for(uword I{1}; I < dt_table.n_rows; ++I)
            if(kappa <= dt_table(I, 0)) {
                out[3] = (dt_table(I, 1) - dt_table(I - 1, 1)) / (dt_table(I, 0) - dt_table(I - 1, 0));
                out[0] = out[3] * (kappa - dt_table(I - 1, 0)) + dt_table(I - 1, 1);
                break;
            }

    const auto backbone = interpolate(kappa, t_table);
    out[4] = backbone[1]; // \md{f}
    out[1] = backbone[0]; // f

    out[2] = out[1] / (1. - out[0]);                                             // \bar{f}
    out[5] = ((1. - out[0]) * out[4] + out[1] * out[3]) * pow(1. - out[0], -2.); // \md{\bar{f}}

    return out;
}

pod6 TableCDP::compute_compression_backbone(const double kappa) const {
    pod6 out;

    if(kappa < dc_table(0, 0)) {
        out[3] = dc_table(0, 1) / dc_table(0, 0); // \md{d}
        out[0] = kappa * out[3];                  // d
    }
    else if(kappa > dc_table(dc_table.n_rows - 1, 0)) {
        out[3] = (1. - dc_table(dc_table.n_rows - 1, 1)) / (1. - dc_table(dc_table.n_rows - 1, 0)); // \md{d}
        out[0] = std::min(1., 1. - (1. - kappa) * out[3]);                                          // d
    }
    else
        for(uword I{1}; I < dc_table.n_rows; ++I)
            if(kappa <= dc_table(I, 0)) {
                out[3] = (dc_table(I, 1) - dc_table(I - 1, 1)) / (dc_table(I, 0) - dc_table(I - 1, 0));
                out[0] = out[3] * (kappa - dc_table(I - 1, 0)) + dc_table(I - 1, 1);
                break;
            }

    const auto backbone = interpolate(kappa, c_table);
    out[4] = backbone[1]; // \md{f}
    out[1] = backbone[0]; // f

    out[2] = out[1] / (1. - out[0]);                                             // \bar{f}
    out[5] = ((1. - out[0]) * out[4] + out[1] * out[3]) * pow(1. - out[0], -2.); // \md{\bar{f}}

    return out;
}

TableCDP::TableCDP(const unsigned T, const double E, const double V, mat&& TT, mat&& CT, mat&& TDT, mat&& CDT, const double GT, const double GC, const double AP, const double BC, const double S, const double R)
    : NonlinearCDP(T, E, V, GT, GC, AP, BC, S, R)
    , t_table(std::move(TT))
    , c_table(std::move(CT))
    , dt_table(std::move(TDT))
    , dc_table(std::move(CDT)) {}

unique_ptr<Material> TableCDP::unique_copy() { return std::make_unique<TableCDP>(*this); }
