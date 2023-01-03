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

#include "Fibre3D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

Fibre3D::Fibre3D(const unsigned T, uvec&& ST)
    : Section3D(T, 0)
    , fibre_tag(std::forward<uvec>(ST)) {}

Fibre3D::Fibre3D(const Fibre3D& old_obj)
    : Section3D(old_obj)
    , fibre_tag(old_obj.fibre_tag) {
    fibre.clear();
    fibre.reserve(old_obj.fibre.size());
    for(const auto& I : old_obj.fibre) fibre.emplace_back(I->get_copy());
}

int Fibre3D::initialize(const shared_ptr<DomainBase>& D) {
    fibre.clear();
    fibre.reserve(fibre_tag.n_elem);

    auto total_a = 0., total_linear_density = 0.;

    for(const auto I : fibre_tag) {
        fibre.emplace_back(suanpan::initialized_section_copy(D, I));
        if(nullptr == fibre.back() || SectionType::D3 != fibre.back()->get_section_type()) {
            SP_W("Section {} is ignored as it is not a valid 3D section.\n", I);
            fibre.pop_back();
        }
        else {
            total_a += fibre.back()->get_parameter(ParameterType::AREA);
            total_linear_density += fibre.back()->get_parameter(ParameterType::LINEARDENSITY);
        }
    }

    access::rw(area) = total_a;
    access::rw(linear_density) = total_linear_density;

    initial_stiffness.zeros();
    access::rw(linear_density) = 0.;
    for(const auto& I : fibre) {
        I->set_eccentricity(I->get_eccentricity() + eccentricity);
        I->initialize_base(D);
        I->initialize(D);
        initial_stiffness += I->get_initial_stiffness();
    }

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Fibre3D::get_copy() { return make_unique<Fibre3D>(*this); }

int Fibre3D::update_trial_status(const vec& t_deformation) {
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

int Fibre3D::clear_status() {
    current_deformation = trial_deformation.zeros();
    current_resistance = trial_resistance.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    auto code = 0;
    for(const auto& I : fibre) code += I->clear_status();
    return code;
}

int Fibre3D::commit_status() {
    current_deformation = trial_deformation;
    current_resistance = trial_resistance;
    current_stiffness = trial_stiffness;
    auto code = 0;
    for(const auto& I : fibre) code += I->commit_status();
    return code;
}

int Fibre3D::reset_status() {
    trial_deformation = current_deformation;
    trial_resistance = current_resistance;
    trial_stiffness = current_stiffness;
    auto code = 0;
    for(const auto& I : fibre) code += I->reset_status();
    return code;
}

void Fibre3D::print() { sp_info("A 3D fibre section that consists of fibres.\n"); }
