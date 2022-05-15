/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "Bar2D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

Bar2D::Bar2D(const unsigned T, const double AR, const unsigned MT, const double EC)
    : Section2D(T, MT, AR, EC) {}

Bar2D::Bar2D(const Bar2D& old_obj)
    : Section2D(old_obj)
    , s_material(suanpan::make_copy(old_obj.s_material)) {}

int Bar2D::initialize(const shared_ptr<DomainBase>& D) {
    s_material = suanpan::initialized_material_copy(D, material_tag);

    access::rw(linear_density) = area * s_material->get_parameter(ParameterType::DENSITY);

    initial_stiffness.set_size(2, 2);
    initial_stiffness(0, 0) = s_material->get_initial_stiffness().at(0) * area;
    initial_stiffness(0, 1) = initial_stiffness(1, 0) = initial_stiffness(0, 0) * eccentricity(0);
    initial_stiffness(1, 1) = initial_stiffness(0, 1) * eccentricity(0);

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Bar2D::get_copy() { return make_unique<Bar2D>(*this); }

int Bar2D::update_trial_status(const vec& t_deformation) {
    trial_deformation = t_deformation;

    if(s_material->update_trial_status(trial_deformation(0) + trial_deformation(1) * eccentricity(0)) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    trial_stiffness(0, 0) = s_material->get_trial_stiffness().at(0) * area;
    trial_stiffness(0, 1) = trial_stiffness(1, 0) = trial_stiffness(0, 0) * eccentricity(0);
    trial_stiffness(1, 1) = trial_stiffness(0, 1) * eccentricity(0);

    trial_resistance(0) = s_material->get_trial_stress().at(0) * area;
    trial_resistance(1) = trial_resistance(0) * eccentricity(0);

    return SUANPAN_SUCCESS;
}

int Bar2D::clear_status() {
    trial_deformation = current_deformation.zeros();
    trial_resistance = current_resistance.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return s_material->clear_status();
}

int Bar2D::commit_status() {
    current_deformation = trial_deformation;
    current_resistance = trial_resistance;
    current_stiffness = trial_stiffness;
    return s_material->commit_status();
}

int Bar2D::reset_status() {
    trial_deformation = current_deformation;
    trial_resistance = current_resistance;
    trial_stiffness = current_stiffness;
    return s_material->reset_status();
}

void Bar2D::print() {
    suanpan_info("A Bar2D section that represents for example rebar in RC section.\n");
    s_material->print();
}
