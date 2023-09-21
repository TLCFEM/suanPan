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

#include "UniversalOS.h"
#include <Domain/DomainBase.h>

UniversalOS::UniversalOS(const unsigned T, const unsigned BT, const unsigned MI, uvec&& FA, uvec&& FB)
    : StressWrapper(T, BT, MI, std::forward<uvec>(FA), std::forward<uvec>(FB), MaterialType::OS) {}

void UniversalOS::print() {
    suanpan_info("An open section material wrapper.\n");
    suanpan_info("Strain:", current_strain);
    suanpan_info("Stress:", current_stress);
    if(base) base->print();
}

OS03::OS03(const unsigned T, const unsigned BT, const unsigned MI)
    : UniversalOS(T, BT, MI, uvec{0, 3}, uvec{1, 2, 4, 5}) {}

unique_ptr<Material> OS03::get_copy() { return make_unique<OS03>(*this); }
