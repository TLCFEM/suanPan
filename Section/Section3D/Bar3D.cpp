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

#include "Bar3D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

Bar3D::Bar3D(const unsigned T, const double AR, const unsigned MT, const double EA, const double EB)
    : Section3D(T, MT, AR, vec{EA, EB}) {}

Bar3D::Bar3D(const Bar3D& old_obj)
    : Section3D(old_obj)
    , s_material(suanpan::make_copy(old_obj.s_material)) {}

int Bar3D::initialize(const shared_ptr<DomainBase>& D) {
    s_material = D->get_material(material_tag)->get_copy();

    access::rw(linear_density) = area * s_material->get_parameter(ParameterType::DENSITY);

    initial_stiffness.set_size(3, 3);

    const auto tmp_a = s_material->get_initial_stiffness().at(0) * area;
    const auto& arm_y = eccentricity(0);
    const auto& arm_z = eccentricity(1);
    initial_stiffness(0, 0) = tmp_a;
    initial_stiffness(1, 1) = tmp_a * arm_y * arm_y;
    initial_stiffness(2, 2) = tmp_a * arm_z * arm_z;
    initial_stiffness(1, 0) = initial_stiffness(0, 1) = tmp_a * arm_y;
    initial_stiffness(2, 0) = initial_stiffness(0, 2) = -tmp_a * arm_z;
    initial_stiffness(2, 1) = initial_stiffness(1, 2) = -tmp_a * arm_y * arm_z;

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Section> Bar3D::get_copy() { return make_unique<Bar3D>(*this); }

int Bar3D::update_trial_status(const vec& t_deformation) {
    trial_deformation = t_deformation;

    const auto& arm_y = eccentricity(0);
    const auto& arm_z = eccentricity(1);

    if(s_material->update_trial_status(trial_deformation(0) + trial_deformation(1) * arm_y - trial_deformation(2) * arm_z) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    const auto tmp_a = s_material->get_trial_stiffness().at(0) * area;
    trial_stiffness(0, 0) = tmp_a;
    trial_stiffness(1, 1) = tmp_a * arm_y * arm_y;
    trial_stiffness(2, 2) = tmp_a * arm_z * arm_z;
    trial_stiffness(1, 0) = trial_stiffness(0, 1) = tmp_a * arm_y;
    trial_stiffness(2, 0) = trial_stiffness(0, 2) = -tmp_a * arm_z;
    trial_stiffness(2, 1) = trial_stiffness(1, 2) = -tmp_a * arm_y * arm_z;
    const auto tmp_b = s_material->get_trial_stress().at(0) * area;
    trial_resistance(0) = tmp_b;
    trial_resistance(1) = tmp_b * arm_y;
    trial_resistance(2) = -tmp_b * arm_z;

    return SUANPAN_SUCCESS;
}

int Bar3D::clear_status() {
    current_deformation.zeros();
    trial_deformation.zeros();
    current_resistance.zeros();
    trial_resistance.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return s_material->clear_status();
}

int Bar3D::commit_status() {
    current_deformation = trial_deformation;
    current_resistance = trial_resistance;
    current_stiffness = trial_stiffness;
    return s_material->commit_status();
}

int Bar3D::reset_status() {
    trial_deformation = current_deformation;
    trial_resistance = current_resistance;
    trial_stiffness = current_stiffness;
    return s_material->reset_status();
}

void Bar3D::print() {
    suanpan_info("A Bar3D section that represents for example rebar in RC section.\n");
    s_material->print();
}
