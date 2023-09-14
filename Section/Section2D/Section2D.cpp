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

#include "Section2D.h"
#include <Material/Material.h>

Section2D::IntegrationPoint::IntegrationPoint(const double C, const double W, unique_ptr<Material>&& M)
    : coor(C)
    , weight(W)
    , s_material(std::forward<unique_ptr<Material>>(M)) {}

Section2D::IntegrationPoint::IntegrationPoint(const IntegrationPoint& old_obj)
    : coor(old_obj.coor)
    , weight(old_obj.weight)
    , s_material(suanpan::make_copy(old_obj.s_material)) {}

Section2D::Section2D(const unsigned T, const unsigned MT, const double A, const double EC)
    : Section(T, SectionType::D2, MT, A, vec{EC, 0.}) {}

int Section2D::update_trial_status(const vec& t_deformation) {
    if(norm((trial_deformation = t_deformation) - current_deformation) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stiffness.zeros();
    trial_resistance.zeros();

    for(const auto& I : int_pt) {
        const auto arm = eccentricity(0) - I.coor;
        if(I.s_material->update_trial_status(trial_deformation(0) + trial_deformation(1) * arm) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        auto tmp_a = I.s_material->get_trial_stiffness().at(0) * I.weight;
        trial_stiffness(0, 0) += tmp_a;
        trial_stiffness(0, 1) += tmp_a *= arm;
        trial_stiffness(1, 1) += tmp_a *= arm;
        tmp_a = I.s_material->get_trial_stress().at(0) * I.weight;
        trial_resistance(0) += tmp_a;
        trial_resistance(1) += tmp_a * arm;
    }

    trial_stiffness(1, 0) = trial_stiffness(0, 1);

    return SUANPAN_SUCCESS;
}

int Section2D::clear_status() {
    current_deformation = trial_deformation.zeros();
    current_resistance = trial_resistance.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->clear_status();
    return code;
}

int Section2D::commit_status() {
    current_deformation = trial_deformation;
    current_resistance = trial_resistance;
    current_stiffness = trial_stiffness;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->commit_status();
    return code;
}

int Section2D::reset_status() {
    trial_deformation = current_deformation;
    trial_resistance = current_resistance;
    trial_stiffness = current_stiffness;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->reset_status();
    return code;
}
