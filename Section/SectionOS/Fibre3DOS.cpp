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

#include "Fibre3DOS.h"
#include <Recorder/OutputType.h>

Fibre3DOS::Fibre3DOS(const unsigned T, uvec&& ST)
    : Fibre(T, std::move(ST), SectionType::OS3D) {}

unique_ptr<Section> Fibre3DOS::get_copy() { return make_unique<Fibre3DOS>(*this); }

std::vector<vec> Fibre3DOS::record(const OutputType P) {
    if(OutputType::S == P) {
        vec force(6, fill::zeros);
        for(const auto& I : fibre) for(const auto& J : I->record(P)) if(!J.empty()) force += J;
        return {force};
    }

    return {};
}
