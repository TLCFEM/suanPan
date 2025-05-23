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

#include "SectionExample.h"

#include <Toolbox/utility.h>

SUANPAN_EXPORT void new_sectionexample(unique_ptr<Section>& return_obj, std::istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    double edge;
    if(!get_input(command, edge)) {
        suanpan_error("A valid width is required.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("A valid modulus is required.\n");
        return;
    }

    return_obj = std::make_unique<SectionExample>(tag, edge, elastic_modulus);
}

SectionExample::SectionExample(const unsigned T, const double S, const double E)
    : Section(T, SectionType::D2, 0, S * S)
    , edge_length(S)
    , moment_inertia(S * S * S * S / 12.)
    , elastic_modulus(E) {}

int SectionExample::initialize(const shared_ptr<DomainBase>&) {
    initial_stiffness.zeros(2, 2);
    initial_stiffness(0, 0) = elastic_modulus * area;
    initial_stiffness(1, 1) = elastic_modulus * moment_inertia;
    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> SectionExample::get_copy() { return std::make_unique<SectionExample>(*this); }

int SectionExample::update_trial_status(const vec& t_deformation) {
    trial_deformation = t_deformation;
    trial_resistance = trial_stiffness * trial_deformation;
    return 0;
}

int SectionExample::clear_status() {
    current_deformation.zeros();
    trial_deformation.zeros();
    current_resistance.zeros();
    trial_resistance.zeros();
    return 0;
}

int SectionExample::commit_status() {
    current_deformation = trial_deformation;
    current_resistance = trial_resistance;
    return 0;
}

int SectionExample::reset_status() {
    trial_deformation = current_deformation;
    trial_resistance = current_resistance;
    return 0;
}

void SectionExample::print() {
    suanpan_info("An example section that represents a square elastic section.\n");
}
