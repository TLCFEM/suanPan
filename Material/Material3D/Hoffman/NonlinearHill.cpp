/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "NonlinearHill.h"

NonlinearHill::NonlinearHill(const unsigned T, vec&& EE, vec&& VV, vec&& S, const double R)
    : NonlinearHoffman(T, std::move(EE), std::move(VV), vec{S(0), S(0), S(1), S(1), S(2), S(2), S(3), S(4), S(5)}, R) {}
