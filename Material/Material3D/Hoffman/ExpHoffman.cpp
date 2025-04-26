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

#include "ExpHoffman.h"

double ExpHoffman::compute_k(const double p_strain) const { return 1. + a - a * exp(-b * p_strain); }

double ExpHoffman::compute_dk(const double p_strain) const { return a * b * exp(-b * p_strain); }

ExpHoffman::ExpHoffman(const unsigned T, vec&& E, vec&& V, vec&& S, const double A, const double B, const double R)
    : DataExpHoffman{A, B}
    , NonlinearHoffman(T, std::move(E), std::move(V), std::move(S), R) {}

unique_ptr<Material> ExpHoffman::get_copy() { return std::make_unique<ExpHoffman>(*this); }
