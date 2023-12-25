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

#include "MaterialTemplate.h"

MaterialTemplate::MaterialTemplate(const unsigned T)
    : Material(T, MaterialType::D1, 0.) {}

int MaterialTemplate::initialize(const shared_ptr<DomainBase>&) { return SUANPAN_SUCCESS; }

double MaterialTemplate::get_parameter(const ParameterType) const { return 0.; }

unique_ptr<Material> MaterialTemplate::get_copy() { return make_unique<MaterialTemplate>(*this); }

int MaterialTemplate::update_trial_status(const vec&) { return 0; }

int MaterialTemplate::clear_status() { return 0; }

int MaterialTemplate::commit_status() { return 0; }

int MaterialTemplate::reset_status() { return 0; }
