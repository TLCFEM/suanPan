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

#include "Fibre1D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

Fibre1D::Fibre1D(const unsigned T, uvec&& ST)
    : Section1D(T, 0, 0.)
    , fibre_tag(std::forward<uvec>(ST)) {}

Fibre1D::Fibre1D(const Fibre1D& old_obj)
    : Section1D(old_obj)
    , fibre_tag(old_obj.fibre_tag) {
    fibre.clear();
    fibre.reserve(old_obj.fibre.size());
    for(const auto& I : old_obj.fibre) fibre.emplace_back(I->get_copy());
}

int Fibre1D::initialize(const shared_ptr<DomainBase>& D) {
    auto total_a = 0., total_linear_density = 0.;

    initial_stiffness.zeros();

    fibre.clear();
    fibre.reserve(fibre_tag.n_elem);

    for(const auto I : fibre_tag) {
        fibre.emplace_back(suanpan::initialized_section_copy(D, I));
        if(nullptr == fibre.back() || SectionType::D1 != fibre.back()->get_section_type()) {
            suanpan_warning("Fibre1D ignores section %llu since it is not a 1D section.\n", I);
            fibre.pop_back();
        }
        else {
            initial_stiffness += fibre.back()->get_initial_stiffness();
            total_a += fibre.back()->get_parameter(ParameterType::AREA);
            total_linear_density += fibre.back()->get_parameter(ParameterType::LINEARDENSITY);
        }
    }

    access::rw(area) = total_a;
    access::rw(linear_density) = total_linear_density;

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Fibre1D::get_copy() { return make_unique<Fibre1D>(*this); }

int Fibre1D::update_trial_status(const vec& t_deformation) {
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

int Fibre1D::clear_status() {
    auto code = 0;
    for(const auto& I : fibre) code += I->clear_status();
    return code;
}

int Fibre1D::commit_status() {
    auto code = 0;
    for(const auto& I : fibre) code += I->commit_status();
    return code;
}

int Fibre1D::reset_status() {
    auto code = 0;
    for(const auto& I : fibre) code += I->reset_status();
    return code;
}

void Fibre1D::print() {
    suanpan_info("A 1D fibre section consists of %llu sections.\n", fibre.size());
    for(const auto& I : fibre) I->print();
}
