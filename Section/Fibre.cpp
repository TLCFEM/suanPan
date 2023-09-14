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

#include "Fibre.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

Fibre::Fibre(const unsigned T, uvec&& FT, const SectionType ST)
    : Section(T, ST, 0)
    , fibre_tag(std::forward<uvec>(FT)) {}

int Fibre::initialize(const shared_ptr<DomainBase>& D) {
    fibre.clear();
    fibre.reserve(fibre_tag.n_elem);

    const auto host_type = get_section_type();

    auto total_area = 0., total_linear_density = 0.;

    initial_stiffness.zeros();
    for(const auto I : fibre_tag) {
        fibre.emplace_back(suanpan::initialized_section_copy(D, I));
        if(nullptr == fibre.back() || host_type != fibre.back()->get_section_type()) {
            suanpan_warning("Section {} is ignored as it is not compatible with fibre section {}.\n", I, get_tag());
            fibre.pop_back();
        }
        else {
            total_area += fibre.back()->get_parameter(ParameterType::AREA);
            total_linear_density += fibre.back()->get_parameter(ParameterType::LINEARDENSITY);
            initial_stiffness += fibre.back()->get_initial_stiffness();
        }
    }

    access::rw(area) = total_area;
    access::rw(linear_density) = total_linear_density;

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int Fibre::update_trial_status(const vec& t_deformation) {
    trial_deformation = t_deformation;

    trial_stiffness.zeros();
    trial_resistance.zeros();

    for(const auto& I : fibre) {
        if(I->update_trial_status(t_deformation) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        trial_stiffness += I->get_trial_stiffness();
        trial_resistance += I->get_trial_resistance();
    }

    return SUANPAN_SUCCESS;
}

int Fibre::clear_status() {
    current_deformation = trial_deformation.zeros();
    current_resistance = trial_resistance.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    auto code = 0;
    for(const auto& I : fibre) code += I->clear_status();
    return code;
}

int Fibre::commit_status() {
    current_deformation = trial_deformation;
    current_resistance = trial_resistance;
    current_stiffness = trial_stiffness;
    auto code = 0;
    for(const auto& I : fibre) code += I->commit_status();
    return code;
}

int Fibre::reset_status() {
    trial_deformation = current_deformation;
    trial_resistance = current_resistance;
    trial_stiffness = current_stiffness;
    auto code = 0;
    for(const auto& I : fibre) code += I->reset_status();
    return code;
}

void Fibre::print() {
    suanpan_info("A composite fibre section consists of the following sections.\n");
    for(const auto& I : fibre) I->print();
}
