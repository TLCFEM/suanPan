/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "Patch.h"
#include <Element/Utility/IGA/BSpline.h>
#include <Material/Material.h>
#include <Section/Section.h>

Patch::Patch(field<vec>&& KT)
    : knot_pool(std::forward<field<vec>>(KT)) {
    element_span.clear();
    element_span.reserve(3);
    for(auto& I : knot_pool) element_span.emplace_back(IGA::compute_all_element_span(I));
}

uvec Patch::get_number_of_control_points() const {
    vector<uword> number;
    for(const auto& I : knot_pool) number.emplace_back(I.n_elem - IGA::compute_order(I) - 1);
    return number;
}

MaterialPatch::MaterialPatch(const unsigned T, const unsigned ND, uvec&& NT, uvec&& MT, field<vec>&& KP, const bool R, const MaterialType MTP)
    : Patch(std::forward<field<vec>>(KP))
    , MaterialElement(T, static_cast<unsigned>(NT.size()), ND, std::forward<uvec>(NT), std::forward<uvec>(MT), R, MTP, {}) {}

MaterialPatch2D::MaterialPatch2D(const unsigned T, const unsigned ND, uvec&& NT, uvec&& MT, field<vec>&& KP, const bool R)
    : MaterialPatch(T, ND, std::forward<uvec>(NT), std::forward<uvec>(MT), std::forward<field<vec>>(KP), R, MaterialType::D2) {}

MaterialPatch3D::MaterialPatch3D(const unsigned T, const unsigned ND, uvec&& NT, uvec&& MT, field<vec>&& KP, const bool R)
    : MaterialPatch(T, ND, std::forward<uvec>(NT), std::forward<uvec>(MT), std::forward<field<vec>>(KP), R, MaterialType::D3) {}

SectionPatch::SectionPatch(const unsigned T, const unsigned ND, uvec&& NT, uvec&& ST, field<vec>&& KP, const bool R, const SectionType STP)
    : Patch(std::forward<field<vec>>(KP))
    , SectionElement(T, static_cast<unsigned>(NT.size()), ND, std::forward<uvec>(NT), std::forward<uvec>(ST), R, STP, {}) {}

SectionPatch2D::SectionPatch2D(const unsigned T, const unsigned ND, uvec&& NT, uvec&& ST, field<vec>&& KP, const bool R)
    : SectionPatch(T, ND, std::forward<uvec>(NT), std::forward<uvec>(ST), std::forward<field<vec>>(KP), R, SectionType::D2) {}

SectionPatch3D::SectionPatch3D(const unsigned T, const unsigned ND, uvec&& NT, uvec&& ST, field<vec>&& KP, const bool R)
    : SectionPatch(T, ND, std::forward<uvec>(NT), std::forward<uvec>(ST), std::forward<field<vec>>(KP), R, SectionType::D3) {}
