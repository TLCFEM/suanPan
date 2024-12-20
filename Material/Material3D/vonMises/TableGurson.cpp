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

#include "TableGurson.h"

vec TableGurson::compute_hardening(const double plastic_strain) const {
    vec response(2);

    for(uword I = 1; I < hardening_table.n_rows; ++I)
        if(hardening_table(I, 0) > plastic_strain) {
            response(1) = (hardening_table(I, 1) - hardening_table(I - 1, 1)) / (hardening_table(I, 0) - hardening_table(I - 1, 0));
            response(0) = hardening_table(I - 1, 1) + (plastic_strain - hardening_table(I - 1, 0)) * response(1);
            return response;
        }

    response(1) = 1E-12;
    response(0) = hardening_table(hardening_table.n_rows - 1, 0);

    return response;
}

TableGurson::TableGurson(const unsigned T, const double E, const double V, mat&& HARDEN, const double Q1, const double Q2, const double FN, const double SN, const double EN, const double R)
    : NonlinearGurson(T, E, V, Q1, Q2, FN, SN, EN, R)
    , hardening_table(std::move(HARDEN)) {}

unique_ptr<Material> TableGurson::get_copy() { return make_unique<TableGurson>(*this); }
