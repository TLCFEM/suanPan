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

#include "MaterialElement.h"
#include <Material/Material.h>

MaterialElement::MaterialElement(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& MT, const bool F, const MaterialType MTP, vector<DOF>&& DI)
    : Element(T, NN, ND, std::move(NT), std::move(MT), F, MTP, std::move(DI)) {}

MaterialElement1D::MaterialElement1D(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& MT, const bool F, vector<DOF>&& DI)
    : MaterialElement(T, NN, ND, std::move(NT), std::move(MT), F, MaterialType::D1, std::move(DI)) {}

MaterialElement2D::MaterialElement2D(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& MT, const bool F, vector<DOF>&& DI)
    : MaterialElement(T, NN, ND, std::move(NT), std::move(MT), F, MaterialType::D2, std::move(DI)) {}

MaterialElement3D::MaterialElement3D(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& MT, const bool F, vector<DOF>&& DI)
    : MaterialElement(T, NN, ND, std::move(NT), std::move(MT), F, MaterialType::D3, std::move(DI)) {}
