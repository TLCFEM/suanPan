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

#include "TableCDP.h"

podarray<double> TableCDP::compute_tension_backbone(const double kappa) const {
    podarray<double> out(6);

    if(kappa < dt_table(0, 0)) {
        out(3) = dt_table(0, 1) / dt_table(0, 0); // \md{d}
        out(0) = kappa * out(3);                  // d
    }
    else if(kappa > dt_table(dt_table.n_rows - 1, 0)) {
        out(3) = (1. - dt_table(dt_table.n_rows - 1, 1)) / (1. - dt_table(dt_table.n_rows - 1, 0)); // \md{d}
        out(0) = std::min(1., 1. - (1. - kappa) * out(3));                                          // d
    }
    else
        for(uword I = 1; I < dt_table.n_rows; ++I) {
            if(kappa <= dt_table(I, 0)) {
                out(3) = (dt_table(I, 1) - dt_table(I - 1, 1)) / (dt_table(I, 0) - dt_table(I - 1, 0));
                out(0) = out(3) * (kappa - dt_table(I - 1, 0)) + dt_table(I - 1, 1);
            }
        }

    out(4) = -t_table(t_table.n_rows - 1, 1) / (1. - t_table(t_table.n_rows - 1, 0)); // \md{f}
    out(1) = std::max(0., (kappa - 1.) * out(4));                                     // f

    for(uword I = 1; I < t_table.n_rows; ++I) {
        if(kappa <= t_table(I, 0)) {
            out(4) = (t_table(I, 1) - t_table(I - 1, 1)) / (t_table(I, 0) - t_table(I - 1, 0));
            out(1) = out(4) * (kappa - t_table(I - 1, 0)) + t_table(I - 1, 1);
        }
    }

    out(2) = out(1) / (1. - out(0));                                             // \bar{f}
    out(5) = ((1. - out(0)) * out(4) + out(1) * out(3)) * pow(1. - out(0), -2.); // \md{\bar{f}}

    return out;
}

podarray<double> TableCDP::compute_compression_backbone(const double kappa) const {
    podarray<double> out(6);

    if(kappa < dc_table(0, 0)) {
        out(3) = dc_table(0, 1) / dc_table(0, 0); // \md{d}
        out(0) = kappa * out(3);                  // d
    }
    else if(kappa > dc_table(dc_table.n_rows - 1, 0)) {
        out(3) = (1. - dc_table(dc_table.n_rows - 1, 1)) / (1. - dc_table(dc_table.n_rows - 1, 0)); // \md{d}
        out(0) = std::min(1., 1. - (1. - kappa) * out(3));                                          // d
    }
    else
        for(uword I = 1; I < dc_table.n_rows; ++I) {
            if(kappa <= dc_table(I, 0)) {
                out(3) = (dc_table(I, 1) - dc_table(I - 1, 1)) / (dc_table(I, 0) - dc_table(I - 1, 0));
                out(0) = out(3) * (kappa - dc_table(I - 1, 0)) + dc_table(I - 1, 1);
            }
        }

    for(uword I = 1; I < dc_table.n_rows; ++I) {
        if(kappa <= dc_table(I, 0)) {
            out(3) = (dc_table(I, 1) - dc_table(I - 1, 1)) / (dc_table(I, 0) - dc_table(I - 1, 0));
            out(0) = out(3) * (kappa - dc_table(I - 1, 0)) + dc_table(I - 1, 1);
        }
    }

    out(4) = -c_table(c_table.n_rows - 1, 1) / (1. - c_table(c_table.n_rows - 1, 0)); // \md{f}
    out(1) = std::min(0., (kappa - 1.) * out(4));                                     // f

    for(uword I = 1; I < c_table.n_rows; ++I) {
        if(kappa <= c_table(I, 0)) {
            out(4) = (c_table(I, 1) - c_table(I - 1, 1)) / (c_table(I, 0) - c_table(I - 1, 0));
            out(1) = out(4) * (kappa - c_table(I - 1, 0)) + c_table(I - 1, 1);
        }
    }

    out(2) = out(1) / (1. - out(0));                                             // \bar{f}
    out(5) = ((1. - out(0)) * out(4) + out(1) * out(3)) * pow(1. - out(0), -2.); // \md{\bar{f}}

    return out;
}

TableCDP::TableCDP(const unsigned T, const double E, const double V, mat&& TT, mat&& CT, mat&& TDT, mat&& CDT, const double AP, const double BC, const double S, const double R)
    : NonlinearCDP(T, E, V, 0., 0., AP, BC, S, R)
    , t_table(std::move(TT))
    , c_table(std::move(CT))
    , dt_table(std::move(TDT))
    , dc_table(std::move(CDT)) {
    t_table.col(0) = abs(t_table.col(0));
    t_table.col(1) = abs(t_table.col(1));
    c_table.col(0) = abs(c_table.col(0));
    c_table.col(1) = -abs(c_table.col(1));
    dt_table = clamp(abs(dt_table), 0., 1.);
    dc_table = clamp(abs(dc_table), 0., 1.);

    const auto TI = t_table.n_rows - 1, TJ = t_table.n_rows - 2;
    const auto CI = c_table.n_rows - 1, CJ = c_table.n_rows - 2;

    const auto lt_slope = (t_table(TI, 1) - t_table(TJ, 1)) / (t_table(TI, 0) - t_table(TJ, 0));
    const auto lc_slope = (c_table(CI, 1) - c_table(CJ, 1)) / (c_table(CI, 0) - c_table(CJ, 0));

    access::rw(g_t) = fabs(trapz(t_table.col(0), t_table.col(1)).eval().at(0)) - t_table(TI, 1) * t_table(TI, 1) / lt_slope;
    access::rw(g_c) = fabs(trapz(c_table.col(0), c_table.col(1)).eval().at(0)) - c_table(TI, 1) * c_table(TI, 1) / lc_slope;

    vec holder = t_table.col(0);
    for(uword I = 0; I < t_table.n_rows; ++I) holder(I) = trapz(t_table.col(0).rows(0, I), t_table.col(1).rows(0, I)).eval().at(0) / g_t;
    t_table.col(0) = holder;

    holder = c_table.col(0);
    for(uword I = 0; I < c_table.n_rows; ++I) holder(I) = trapz(c_table.col(0).rows(0, I), c_table.col(1).rows(0, I)).eval().at(0) / g_c;
    c_table.col(0) = holder;
}

unique_ptr<Material> TableCDP::get_copy() { return make_unique<TableCDP>(*this); }
