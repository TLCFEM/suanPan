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

#include "Trivial.h"

Trivial::Trivial(const unsigned T)
    : Material1D(T, 0.) {}

int Trivial::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = 0.;
    trial_stress = current_stress = 0.;

    ConstantStiffness(this);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Trivial::get_copy() { return make_unique<Trivial>(*this); }

int Trivial::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    return SUANPAN_SUCCESS;
}

int Trivial::clear_status() {
    current_strain.zeros();
    return reset_status();
}

int Trivial::commit_status() {
    current_strain = trial_strain;
    return SUANPAN_SUCCESS;
}

int Trivial::reset_status() {
    trial_strain = current_strain;
    return SUANPAN_SUCCESS;
}

void Trivial::print() {
    suanpan_info("A trivial material model returns zero response.\n");
}
