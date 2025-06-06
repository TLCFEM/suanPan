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

#include "Combine.h"

#include <Domain/DomainBase.h>

Combine::Combine(const unsigned T, uvec&& TP)
    : Amplitude(T)
    , tag_pool(std::move(TP)) {}

void Combine::initialize(const shared_ptr<DomainBase>& D) {
    for(const auto I : tag_pool)
        if(D->find<Amplitude>(I)) amp_pool.emplace_back(D->initialized_amplitude_copy(I));
}

unique_ptr<Amplitude> Combine::get_copy() { return std::make_unique<Combine>(*this); }

double Combine::get_amplitude(const double T) {
    auto A = 1.;
    for(const auto& amp : amp_pool) A *= amp->get_amplitude(T);
    return A;
}

void Combine::print() {
    suanpan_info("Combine.\n");
}
