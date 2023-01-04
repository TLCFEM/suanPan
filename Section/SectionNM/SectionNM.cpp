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

#include "SectionNM.h"

void SectionNM::initialize_history(const unsigned size) {
    if(initial_history.empty()) initial_history.zeros(size);
    else if(static_cast<uword>(size) > initial_history.size()) initial_history.resize(size);

    trial_history = current_history = initial_history;
}

int SectionNM::clear_status() {
    current_deformation = trial_deformation.zeros();
    current_resistance = trial_resistance.zeros();
    current_history = trial_history.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    return SUANPAN_SUCCESS;
}

int SectionNM::commit_status() {
    current_deformation = trial_deformation;
    current_resistance = trial_resistance;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int SectionNM::reset_status() {
    trial_deformation = current_deformation;
    trial_resistance = current_resistance;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void SectionNM::print() {
    suanpan_info("An N-M interaction based section.\n");
}
