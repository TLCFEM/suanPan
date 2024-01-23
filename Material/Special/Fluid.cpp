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

#include "Fluid.h"

Fluid::Fluid(const unsigned T, const double E, const double R)
    : DataFluid{E}
    , Material(T, MaterialType::DS, R) {}

int Fluid::initialize(const shared_ptr<DomainBase>&) { return SUANPAN_SUCCESS; }

double Fluid::get_parameter(const ParameterType P) const {
    if(ParameterType::BULKMODULUS == P) return bulk_modulus;
    return 0.;
}

unique_ptr<Material> Fluid::get_copy() { return make_unique<Fluid>(*this); }

int Fluid::update_trial_status(const vec&) { return SUANPAN_SUCCESS; }

int Fluid::clear_status() { return SUANPAN_SUCCESS; }

int Fluid::commit_status() { return SUANPAN_SUCCESS; }

int Fluid::reset_status() { return SUANPAN_SUCCESS; }
