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

#include "SectionElement.h"
#include <Section/Section.h>

SectionElement::SectionElement(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& ST, const bool F, const SectionType STP)
    : Element(T, NN, ND, std::forward<uvec>(NT), std::forward<uvec>(ST), F, STP) {}

SectionElement1D::SectionElement1D(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& ST, const bool F)
    : SectionElement(T, NN, ND, std::forward<uvec>(NT), std::forward<uvec>(ST), F, SectionType::D1) {}

SectionElement2D::SectionElement2D(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& ST, const bool F)
    : SectionElement(T, NN, ND, std::forward<uvec>(NT), std::forward<uvec>(ST), F, SectionType::D2) {}

SectionElement3D::SectionElement3D(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& ST, const bool F)
    : SectionElement(T, NN, ND, std::forward<uvec>(NT), std::forward<uvec>(ST), F, SectionType::D3) {}

SectionNMElement2D::SectionNMElement2D(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& ST, const bool F)
    : SectionElement(T, NN, ND, std::forward<uvec>(NT), std::forward<uvec>(ST), F, SectionType::NM2D) {}

SectionNMElement3D::SectionNMElement3D(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& ST, const bool F)
    : SectionElement(T, NN, ND, std::forward<uvec>(NT), std::forward<uvec>(ST), F, SectionType::NM3D) {}
