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
            }

    out[4] = -t_table(t_table.n_rows - 1, 1) / (1. - t_table(t_table.n_rows - 1, 0)); // \md{f}
    out[1] = std::max(0., (kappa - 1.) * out[4]);                                     // f

    for(uword I{1}; I < t_table.n_rows; ++I)
        if(kappa <= t_table(I, 0)) {
            out[4] = (t_table(I, 1) - t_table(I - 1, 1)) / (t_table(I, 0) - t_table(I - 1, 0));
            out[1] = out[4] * (kappa - t_table(I - 1, 0)) + t_table(I - 1, 1);
        }

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
            }

    for(uword I{1}; I < dc_table.n_rows; ++I)
        if(kappa <= dc_table(I, 0)) {
            out[3] = (dc_table(I, 1) - dc_table(I - 1, 1)) / (dc_table(I, 0) - dc_table(I - 1, 0));
            out[0] = out[3] * (kappa - dc_table(I - 1, 0)) + dc_table(I - 1, 1);
        }

    out[4] = -c_table(c_table.n_rows - 1, 1) / (1. - c_table(c_table.n_rows - 1, 0)); // \md{f}
    out[1] = std::min(0., (kappa - 1.) * out[4]);                                     // f

    for(uword I{1}; I < c_table.n_rows; ++I)
        if(kappa <= c_table(I, 0)) {
            out[4] = (c_table(I, 1) - c_table(I - 1, 1)) / (c_table(I, 0) - c_table(I - 1, 0));
            out[1] = out[4] * (kappa - c_table(I - 1, 0)) + c_table(I - 1, 1);
        }

    out[2] = out[1] / (1. - out[0]);                                             // \bar{f}
    out[5] = ((1. - out[0]) * out[4] + out[1] * out[3]) * pow(1. - out[0], -2.); // \md{\bar{f}}

    return out;
}

TableCDP::TableCDP(const unsigned T, const double E, const double V, mat&& TT, mat&& CT, mat&& TDT, mat&& CDT, const double AP, const double BC, const double S, const double R)
    : NonlinearCDP(T, E, V, 0., 0., AP, BC, S, R)
    , t_table(abs(std::move(TT)))
    , c_table(abs(std::move(CT)))
    , dt_table(abs(std::move(TDT)))
    , dc_table(abs(std::move(CDT))) {
    const auto convert = [](mat& backbone, mat& damage, const double g) {
        // convert plastic strain (first column) to kappa
        vec kappa = backbone.col(0);
        for(uword I{0}; I < backbone.n_rows; ++I) kappa(I) = as_scalar(trapz(backbone.col(0).rows(0, I), backbone.col(1).rows(0, I))) / g;

        // convert damage table
        vec d_kappa;
        interp1(backbone.col(0), kappa, damage.col(0), d_kappa);
        damage.col(0) = d_kappa;
        damage.col(1).clamp(0., 1.);

        // check-in kappa for backbone table
        backbone.col(0) = kappa;
    };

    convert(t_table, dt_table, access::rw(g_t) = as_scalar(trapz(t_table.col(0), t_table.col(1))));
    convert(c_table, dc_table, access::rw(g_c) = as_scalar(trapz(c_table.col(0), c_table.col(1))));

    c_table.col(1) *= -1.; // ensure compression table is negative
}

unique_ptr<Material> TableCDP::unique_copy() { return std::make_unique<TableCDP>(*this); }
