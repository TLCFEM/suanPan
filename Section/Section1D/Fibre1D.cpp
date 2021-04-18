/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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

Fibre1D::Fibre1D(const unsigned T, podarray<unsigned>&& ST)
	: Section1D(T, 0, 0.)
	, fibre_tag(std::forward<podarray<unsigned>>(ST)) {}

Fibre1D::Fibre1D(const Fibre1D& old_obj)
	: Section1D(old_obj)
	, fibre_tag(old_obj.fibre_tag) {
	fibre.clear();
	fibre.reserve(old_obj.fibre.size());
	for(const auto& I : old_obj.fibre) fibre.emplace_back(I->get_copy());
}

void Fibre1D::initialize(const shared_ptr<DomainBase>& D) {
	auto total_a = 0.;
	vec total_q(2, fill::zeros);

	fibre.clear();
	fibre.reserve(fibre_tag.n_elem);

	for(uword I = 0; I < fibre_tag.n_elem; ++I)
		if(D->find_section(fibre_tag(I))) {
			if(const auto& t_section = D->get_section(fibre_tag(I)); SectionType::D1 == t_section->get_section_type()) {
				if(!t_section->is_initialized()) {
					t_section->Section::initialize(D);
					t_section->initialize(D);
					t_section->set_initialized(true);
				}
				fibre.emplace_back(t_section->get_copy());
				const auto t_area = fibre.back()->get_parameter(ParameterType::AREA);
				total_a += t_area;
				total_q += fibre.back()->get_eccentricity() * t_area;
			}
			else suanpan_warning("skip Section %u since it is not a 1D section.\n", fibre_tag(I));
		}

	access::rw(area) = total_a;
	access::rw(eccentricity) = total_q / total_a;

	initial_stiffness.zeros();
	linear_density = 0.;
	for(const auto& I : fibre) {
		I->set_eccentricity(I->get_eccentricity() - eccentricity);
		I->Section::initialize(D);
		I->initialize(D);
		initial_stiffness += I->get_initial_stiffness();
		linear_density += I->get_parameter(ParameterType::LINEARDENSITY);
	}

	trial_stiffness = current_stiffness = initial_stiffness;
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
	suanpan_info("A 1D fibre section consists of %lu sections.\n", fibre.size());
	for(const auto& I : fibre) I->print();
}
